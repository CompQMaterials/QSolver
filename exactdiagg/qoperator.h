#ifndef QOPERATOR_H
#define QOPERATOR_H

#include <vector>
#include <cmath>
#include <functional>
//#include<omp.h>
#include "fockbasis.h"

//-----------------------------------------------------------------------FermiOp---------------------------

struct FermiOp
{
    int pos;
    bool dag;

    FermiOp(int pos,bool dag):pos(pos),dag(dag){}

    //template<int L>
    bool GiveNull(const FockState& f) const {return f.test(pos)==dag;}

    //template<int L>
    bool ApplyTo(FockState& f) const
    {
        f.flip(pos);
        return BitParityLeft(f,pos);
    }
};

inline FermiOp Create(int pos) {return FermiOp(pos,true);}
inline FermiOp Destroy(int pos) {return FermiOp(pos,false);}


//--------------------------------------------------------------------- StringOp ----------------------------------


class StringOp
{
    std::vector<FermiOp> ops;
public:
    StringOp(){}
    StringOp(const FermiOp& o) {ops.push_back(o);}

    void operator*=(const FermiOp& o) {ops.push_back(o);}
    StringOp operator*(const FermiOp& o) const
    {
        StringOp op=*this;
        op*=o;
        return op;
    }

    void operator *=(const StringOp& so)
    {
       for(uint i = 0; i < so.ops.size(); i++)
          ops.push_back(so.ops[i]);
    }
    StringOp operator *(const StringOp& so)
    {
        StringOp op=*this;
        op *= so;
        return op;
    }
    //template<int L>
    int ApplyTo(FockState& f) const
    {
        int sg=1;
        for(int i=ops.size()-1; i>=0; i--)
        {
            if (ops[i].GiveNull(f))
                return 0;
            bool parity=ops[i].ApplyTo(f);
            if (parity) sg=-sg;
        }
        return sg;
    }
};

inline StringOp operator*(const FermiOp& o1,const FermiOp& o2)
{
    StringOp op=o1;
    op*=o2;
    return op;
}


//--------------------------------------------------------------------- QOperator ----------------------------------

//#include "binder.h"
#include "symmetrygroup.h"

template<class T=double>
struct QOperatorG
{
    std::vector<StringOp> ops;
    std::vector<T> coeff;

    QOperatorG(){}
    QOperatorG(const StringOp& op) {ops.push_back(op); coeff.push_back(1);}

    void Add(const StringOp& op, T c)
    {
        ops.push_back(op);
        coeff.push_back(c);
    }
    void operator+=(const QOperatorG& B)
    {
        for(uint i=0;i<B.coeff.size();i++)
            Add(B.ops[i],B.coeff[i]);
    }
    QOperatorG operator+(const QOperatorG& B) {QOperatorG r=*this;r+=B;return r;}

    void operator*=(T c) {for(T& x:coeff) x*=c;}
    QOperatorG operator*(T c) const {QOperatorG r=*this;r*=c;return r;}

    template<class scalar>
    State<scalar> ApplyTo(const FockState& f) const
    {
        State<scalar> state;
        for(uint o=0;o<ops.size();o++)
        {
            FockState c1=f;
            int sg=ops[o].ApplyTo(c1);
            if (!sg) continue;
            if (sg==1) state[c1]+=coeff[o];
            else       state[c1]-=coeff[o];
        }
        return state;
    }

    template<class scalar>
    State<scalar> ApplyTo(const State<scalar>& s) const
    {
        State<scalar> state;
        for(auto it:s)
        {
            State<scalar> tmp=ApplyTo<scalar>(it.first);
            for(auto it2:tmp)
                state[it2.first]+=it.second*it2.second;
        }
        return state;
    }

    template<class scalar=double>
    SpMat<scalar> toMatrix(const FockBasis& b) const
    {
        SpMat<scalar> M(b.Size(),b.Size());
        int n = 0;
        for(int i=0;i<b.Size();i++)
        {
            State<scalar> res = ApplyTo<scalar>(b.vec[i]);

            for(const auto& it:res)
            {
                int pos=b.PosOf(it.first);
                if(pos  == -1)
                {
                    //cout<<"State "<<it.first<<" is not in the basis"<<endl;
                    n++;
                    //std::terminate();
                }
                else //if( pos ==0 || pos==i || i==0)
                    M(pos,i)=it.second;
            }

        }
        //cout<<n<<" missing states in the basis"<<endl;
        return M;
    }

    template<class scalar=double>
    SpMat<scalar> toMatrix(FockBasis &b, bool add = false, double tol = 1e-10) const
    {
        SpMat<scalar> M(b.Size(),b.Size());
        int n = 0;
        for(int i = 0;i < b.Size(); i++)
        {
            State<scalar> res=ApplyTo<scalar>(b.vec[i]);

            for(const auto& it:res)
            {
                if(fabs(it.second) > tol)
                {
                int pos=b.PosOf(it.first);
                if(pos  == -1)
                {
                 //   cout<<"State "<<it.first<<" is not in the basis"<<endl;
                    n++;
//                    if(add)
//                    b.Add(it.first, 1.0);
//                    if(n>10)
//                    std::terminate();
//                    cout<<fabs(it.second)<<endl;
                }
                else M(pos,i)=it.second;
                }
            }

        }
        if(n>0)
        cout<<n<<" :missing states in the basis"<<endl;
        return M;
    }

    template<class scalar=double>
    SpMat<scalar> toMatrix(const FockBasis &b1,const FockBasis &b2) const
    {
        SpMat<scalar> M(b1.Size(),b2.Size());
        for(int i = 0; i < b2.Size(); i++)
        {
            State<scalar> res=ApplyTo<scalar>(b2.vec[i]);
            for(const auto& it:res)
            {
                int pos=b1.PosOf(it.first);
                if(pos  != -1)
                M(pos,i)=it.second;
            }
        }
        return M;
    }

    template<class scalar>
    SpMat<scalar> toMatrix(const FockBasis &b, const SymmetryGroup<scalar>& G) const
    {

        SpMat<scalar> M(b.Size(),b.Size());
        for(int i=0;i<b.Size();i++)
        {
            if (b.norma[i]==0.0) continue;
            State<scalar> res=ApplyTo<scalar>(b.vec[i]);
            for(const auto& it:res)
            {
                FockState r=it.first;
                scalar gm= G.Representant(r,b.sym) ;
                int pos=b.PosOf(r);
                if (pos!=-1 && b.norma[pos]!=0.0)
                    M(pos,i)+=it.second*gm*b.norma[pos]*b.norma[i];
            }
        }
        return M;
    }

    template<class scalar=double>
    SpMat<scalar> toMatrix(const FockBasis &b1, const SymmetryGroup<scalar>& G,const FockBasis &b2) const
    {
        SpMat<scalar> M(b1.Size(),b2.Size());
        for(int i=0;i<b2.Size();i++)
        {
            if  (b2.norma[i]==0.0) continue;
            State<scalar> res=ApplyTo<scalar>(G.ApplyTo(b2.vec[i],b2.sym));
            for(const auto& it:res)
            {
                FockState r=it.first;
                scalar gm=G.Representant(r,b1.sym);
                int pos=b1.PosOf(r);
                if (pos!=-1 && b1.norma[pos]!=0.0)
                    M(pos,i)+=gm*it.second*b1.norma[pos]*b2.norma[i];
            }
        }
        return M;
    }
};

using QOperator = QOperatorG<double>;
using QOperatorC = QOperatorG< std::complex<double> >;

#endif // QOPERATOR_H
