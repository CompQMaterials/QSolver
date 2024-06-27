#ifndef SYMMETRYGROUP_H
#define SYMMETRYGROUP_H
#include<functional>
#include <armadillo>
#include <complex>

#include "fockbasis.h"
//#include "fockstate.h"





//---------------------------------------------------------- ElementaryOp -------------------------------------




using cmpx=std::complex<double>;

using namespace arma;


using ElementaryOp = std::function<int(FockState&)> ;


inline int IdentityOp(FockState& f){return 1;}


struct TranslationOp
{
    int delta;

    TranslationOp(int delta) : delta(delta) {}

    int operator()(FockState& f) const
    {
        int L = f.size();
        if (delta == 0) return 1;
        FockState l = f >> delta;
        FockState r = f << (L - delta);
        f = l | r;
        bool parity = (r.count() * l.count()) & 1;
        return parity ? -1 : 1;
    }

    TranslationOp Pow(int n) const { return TranslationOp(delta * n); }
};


inline int ReflectionOp(FockState& f)
{
    std::string s;
    to_string(f, s);
    std::reverse(s.begin(), s.end());
    f = FockState(s);
    return (f.count() / 2) & 1 ? -1 : 1;
}


inline int ParticleHoleOp(FockState& f)
{
    f.flip();
    return f.count() & 1 ? -1 : 1;
}


struct ProductOp
{
    std::function<int(FockState&)> op1;
    std::function<int(FockState&)> op2;

    ProductOp(std::function<int(FockState&)> op1, std::function<int(FockState&)> op2) : op1(op1), op2(op2) {}

    int operator()(FockState& f) const
    {
        return (op2(f) == -1) ? -op1(f) : op1(f);
    }
};

template<int R>
struct TensorPowOp
{
    std::function<int(FockState&)> op;
    TensorPowOp(std::function<int(FockState&)> op) : op(op) {}

    int operator()(FockState& f) const
    {
        int L = f.size();
        auto a = Split(f,R);
        int sg = 1;
        for (int i = 0; i < R; i++)
        {
            if (op(a[i]) == -1) sg = -sg;
        }
        f = Join(a);
        return sg;
    }

    TensorPowOp Pow(int n) const { return TensorPowOp(op); }
};

template<int R>
TensorPowOp<R> TensorPow(std::function<int(FockState&)> op) { return TensorPowOp<R>(op); }

template<class Oper, class scalar = double>
std::map<FockState, scalar> ApplyOp(const Oper& Op, const std::map<FockState, scalar>& s)
{
    std::map<FockState, scalar> sr;
    for (const auto& it : s)
    {
        FockState r = it.first;
        scalar coef = it.second;
        int sg = Op(r);
        sr[r] += coef * scalar(sg);
    }
    return sr;
}

//---------------------------------------------------------- SymmetryGroup ----------------------------



template<class scalar>
class SymmetryGroup
{
public:
    std::vector<ElementaryOp> T;
    std::vector<Col<scalar>> Gamma;

    void Add(ElementaryOp op, const Col<scalar>& coeff)
    {
        T.push_back(op);
        Gamma.push_back(coeff);
    }

    int nSym() const { return Gamma.empty() ? 0 : Gamma[0].size(); }

    State<scalar> ApplyTo(const FockState & f, int sym) const
    {
        std::vector<FockState> vec(T.size());
        std::vector<int> sg(T.size());

        for (size_t g = 0; g < T.size(); ++g)
        {
            vec[g] = f;
            sg[g] = T[g](vec[g]);
        }

        State<scalar> state;

        for (size_t g = 0; g < T.size(); ++g)
        {
            scalar& value = state[vec[g]];

            if (sg[g] == 1)
                value += Gamma[g][sym];
            else
                value -= Gamma[g][sym];
        }

        return state;
    }

    State<scalar> ApplyTo(const State<scalar>& s, int sym) const
    {
        State<scalar> state;

        for (const auto& it : s)
        {
            State<scalar> tmp = ApplyTo(it.first, sym);

            for (const auto& it2 : tmp)
                state[it2.first] += it.second * it2.second;
        }

        return state;
    }

    scalar Representant(FockState & f, int sym) const
    {
        State<scalar> res = ApplyTo(f, sym);

        if (res.empty())
        {
            std::cerr << "Representant not found!";
            return 0;
        }

        f = res.begin()->first;
        return res.begin()->second;
    }

    template<class scalar2>
    SymmetryGroup<scalar> DirectProd(const SymmetryGroup<scalar2>& G2) const
    {
        const SymmetryGroup& G1 = *this;

        if (G1.nSym() == 0)
            throw std::invalid_argument("DirectProd: invalid G1");

        SymmetryGroup<scalar> G;

        for (size_t i = 0; i < G1.Gamma.size(); ++i)
        {
            for (size_t j = 0; j < G2.Gamma.size(); ++j)
            {
                auto g = [&](FockState& f) {
                    return (G2.T[j](f) == -1) ? -G1.T[i](f) : G1.T[i](f);
                };

                G.Add(g, kron(G1.Gamma[i], G2.Gamma[j]));
            }
        }

        return G;
    }
};

template<class ElementaryOp>
SymmetryGroup<cmpx> CyclicGroupPow(const ElementaryOp& T, int nElem)
{
    SymmetryGroup<cmpx> G;

    for (int i = 0; i < nElem; ++i)
    {
        Col<cmpx> coeff(nElem);

        for (int k = 0; k < nElem; ++k)
        {
            cmpx phase = {0.0, -2 * M_PI * i * k / nElem};
            coeff[k] = std::exp(phase) * (1.0 / nElem);
        }

        G.Add(T.Pow(i), coeff);
    }

    return G;
}

template<class scalar = double, class ElementaryOp>
SymmetryGroup<scalar> Z2_Group(const ElementaryOp& T)
{
    SymmetryGroup<scalar> G;
    G.Add(IdentityOp, {0.5, 0.5});
    G.Add(T, {0.5, -0.5});
    return G;
}

template<class Oper, class scalar = double>
SpMat<scalar> toMatrix(const Oper& Op, const FockBasis& b1, const SymmetryGroup<scalar>& G, const FockBasis& b2)
{
    SpMat<scalar> M(b1.Size(), b2.Size());
    for (int i = 0; i < b2.Size(); i++)
    {
        if (b2.norma[i] == 0.0) continue;
        auto ssym = G.ApplyTo(b2.vec[i], b2.sym);
        auto res = ApplyOp<Oper, scalar>(Op, ssym);
        for (const auto& it : res)
        {
            FockState r = it.first;
            scalar gm = G.Representant(r, b1.sym);
            int pos = b1.PosOf(r);
            if (pos != -1 && b1.norma[pos] != 0.0)
                M(pos, i) += gm * it.second * b1.norma[pos] * b2.norma[i];
        }
    }
    return M;
}





//template<int L, class scalar>
//class SymmetryGroup
//{
//public:
//    std::vector< ElementaryOp<L> > T;
//    std::vector< Col<scalar> > Gamma;

//    void Add(ElementaryOp<L> op, const Col<scalar>& coeff)
//    {
//        T.push_back(op);
//        Gamma.push_back(coeff);
//    }

//    int nSym() const {return Gamma.size()==0?0:Gamma[0].size();}

//    State<L,scalar> ApplyTo(const FockState& f,int sym) const
//    {
//        std::vector<FockState> vec(T.size());
//        std::vector<int> sg(T.size());
////        #pragma omp parallel for
//        for(uint g=0; g<T.size(); g++)
//        {
//            vec[g]=f;
//            sg[g]=T[g](vec[g]);
//        }
//        State<L,scalar> state;
//        for(uint g=0;g<T.size();g++)
//        {
//            scalar& value=state[vec[g]];
//            if (sg[g]==1) value+=Gamma[g][sym];
//            else          value-=Gamma[g][sym];
////            if (value==0.0) state.erase(vec[g]);
//        }
//        return state;
//    }

//    State<L,scalar> ApplyTo(const State<L,scalar>& s,int sym) const
//    {
//        State<L,scalar> state;
//        for(auto it:s)
//        {
//            State<L,scalar> tmp=ApplyTo(it.first,sym);
//            for(auto it2:tmp)
//                state[it2.first]+=it.second*it2.second;
//        }
//        return state;
//    }

//    scalar Representant(FockState<L>& f, int sym) const
//    {
//        State<L,scalar> res=ApplyTo(f,sym);
//        if (res.empty()) {std::cerr<<"Representant not found!"; return 0;}
//        f=res.begin()->first;
//        return res.begin()->second;
//    }

//    template<class scalar2>
//    SymmetryGroup DirectProd(const SymmetryGroup<L,scalar2>& G2) const
//    {
//        const SymmetryGroup& G1=*this;
//        if (G1.nSym()==0)
//            throw std::invalid_argument("DirectProd: invalid G1");
//        SymmetryGroup G;
//        for(uint i=0;i<G1.Gamma.size();i++)
//            for(uint j=0;j<G2.Gamma.size();j++)
//            {
//                auto g = [=](FockState<L>& f) { return (G2.T[j](f)==-1)?-G1.T[i](f):G1.T[i](f); };
////                auto g=ProductOp<L>(G1.T[i],G2.T[j]);
//                G.Add(g,kron(G1.Gamma[i],G2.Gamma[j]));
//            }
////         G.SetSym( G1.currSym*G2.nSym()+G2.currSym );
//         return G;
//    }
//};

//template<int L,class ElementaryOp>
//SymmetryGroup<L,cmpx> CyclicGroupPow(const ElementaryOp& T,int nElem)
//{
//    SymmetryGroup<L,cmpx> G;
//    for(int i=0;i<nElem;i++)
//    {
//        Col<cmpx> coeff(nElem);
//        for(int k=0;k<nElem;k++)
//        {
//            cmpx phase={0.0, -2*M_PI*i*k/nElem};
//            coeff[k]=std::exp(phase) * (1.0/nElem) ;
//        }
//        G.Add(T.Pow(i),coeff);
//    }
//    return G;
//}

////#include "fockbasis.h"

//template<int L, class scalar=double,class ElementaryOp>
//SymmetryGroup<L,scalar> Z2_Group(const ElementaryOp& T)
//{
//    SymmetryGroup<L,scalar> G;
//    G.Add(IdentityOp<L>,{0.5,0.5});
//    G.Add(T,{0.5,-0.5});
//    return G;
//}



#endif // SYMMETRYGROUP_H
