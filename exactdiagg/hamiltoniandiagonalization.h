#ifndef HAMILTONIANDIAGONALIZATION_H
#define HAMILTONIANDIAGONALIZATION_H
//#include <iostream>
//#include <complex.h>
//#include <algorithm>
//#include <valarray>
//#include <cmath>
//#include<functional>

#include "qoperator.h"
#include "hamsolver.h"

using namespace std;





inline EigenState FindGS(QOperator& ham, FockBasis& b,int nPart, bool add = false)
{
    auto H=ham.toMatrix(b,add);
    cout<<"matrix nonzero="<<H.n_nonzero<<endl; cout.flush();

    vec eval;
    mat evec;
    if(b.Size()>10)
        eigs_sym(eval,evec,H ,std::min(b.Size()-2,101),"sa");
    else
        eig_sym(eval,evec,mat(H));
    uvec ind = sort_index(eval);

    //    for(uint i=0;i<eval.size();i++)
    //        if (eval(i)-eval(ind(0))<1e-10)
    //            cout<<eval[i]<<"\n";
    std::vector<double> unique_eval={ eval(ind(0)) };
    for(uint i=0;i<eval.size();i++)
    {
        const auto x=eval(ind(i));
        if (x-unique_eval.back()>1e-10)
            unique_eval.push_back(x);
        if (unique_eval.size()>2) break; // only print the first 2 different evals
        cout<<x<<"\n";
    }
    cout<<"dim="<<b.Size()<<" nPart="<<nPart<<" sym="<<b.sym<<" ener="<<eval(ind(0))<<endl;
    cout<<"norm:"<<norm(evec.col(ind(0)))<<endl;
    return {eval(ind(0)),evec.col(ind(0)),nPart,b.sym};

    //    Col<double> wf(b.Size());
    //    wf.randu();
    //    auto sol=Diagonalize<double>(H,wf);
    //    cout<<sol.cIter<<" lanczos iter; "<<"nPart="<<nPart<<" sym="<<b.sym<<" ener="<<sol.ener0<<endl;
    //    return {sol.ener0,sol.x0,nPart,b.sym};
}


inline vector<EigenState> FindAllStates(QOperator& ham, FockBasis& b,int nPart, int n_states = 100)
{
    auto H = ham.toMatrix(b);
    //cout<<"matrix nonzero="<<H.n_nonzero<<endl; cout.flush();
    cout<<"matrix nonzero="<<H.n_nonzero<<" Hermitian(tol=1e-6): "<<H.is_hermitian(1e-6)<<endl; cout.flush();
    //mat(H).print("");
    vec eval;
    mat evec;
    if(b.Size() > 5000)
        //eigs_sym(eval, evec, H , 1, "sa");
        eigs_sym(eval,evec,H ,std::min(b.Size()-1,n_states),"sa");
    else
        eig_sym(eval,evec,mat(H));
    uvec ind = sort_index(eval);
    std::vector<double> unique_eval={ eval(ind(0)) };
    vector<EigenState> res;
    for(uint i=0;i<eval.size();i++)
    {
        // cout<<eval(i)<<endl;
        res.push_back(EigenState {eval(ind(i)),evec.col(ind(i)),nPart,b.sym} );
    }

    return  res;
}
using EigenStateC = EigenStateG< std::complex<double> >;

inline vector<EigenStateC> FindAllStates(QOperatorC& ham, FockBasis& b, int nPart, int n_states, double tol = 1e-30, int DenseLim = 1000)
{
    const bool add = false;
    auto H = ham.toMatrix<std::complex<double> >(b, add, tol);
    //H = 0.5*(H + H.t());
    cout<<"matrix nonzero="<<H.n_nonzero<<" Hermitian(tol=1e-4): "<<H.is_hermitian(1e-4)<<endl; cout.flush();
    //H.print("");
    cx_vec eval;
    cx_mat evec;
    if(b.Size() > DenseLim)
        //eigs_sym(eval, evec, H , 1, "sa");
        eigs_gen(eval,evec, H, std::min(b.Size()-2,n_states),"sr");
    else
        eig_gen(eval,evec,cx_mat(H));
    uvec ind = sort_index(real(eval));
    //std::vector< std::complex<double> > unique_eval={ eval(ind(0)) };
    vector<EigenStateC> res;
    for(uint i = 0; i < eval.size(); i++)
    {
        // cout<<eval(i)<<endl;
        res.push_back(EigenStateC {real(eval(ind(i))),evec.col(ind(i)),nPart,b.sym} );
        //res.push_back(EigenStateC {real(eval(i)),evec.col(i),nPart,b.sym} );
    }

    return  res;
}

using EigenStateC = EigenStateG< std::complex<double> >;

inline vector<EigenStateC> FindAllStatesDense(QOperatorC& ham, FockBasis& b, int nPart, vec &eval, cx_mat &evec, double tol = 1e-30)
{
    const bool add = false;
    auto H = ham.toMatrix<std::complex<double> >(b, add, tol);
    //H = 0.5*(H + H.t());
    cout<<"matrix nonzero="<<H.n_nonzero<<" Hermitian(tol=1e-4): "<<H.is_hermitian(1e-4)<<endl; cout.flush();
    //H.print("");


    eig_sym(eval,evec,cx_mat(H));
    uvec ind = sort_index(real(eval));
    //std::vector< std::complex<double> > unique_eval={ eval(ind(0)) };
    vector<EigenStateC> res;
    for(uint i = 0; i < eval.size(); i++)
        res.push_back(EigenStateC {real(eval(ind(i))),evec.col(ind(i)),nPart,b.sym} );

    return  res;
}


//using EigenStateC = EigenStateG< std::complex<double> >;

inline vector<EigenStateC> FindAllStatesDense(QOperatorC& ham, FockBasis& b, int nPart, double tol = 1e-30)
{
    const bool add = false;
    auto H = ham.toMatrix<std::complex<double> >(b, add, tol);
    //H = 0.5*(H + H.t());
    cout<<"matrix nonzero="<<H.n_nonzero<<" Hermitian(tol=1e-4): "<<H.is_hermitian(1e-4)<<endl; cout.flush();
    //H.print("");

    vec eval;
    cx_mat evec;
    eig_sym(eval,evec,cx_mat(H));
    uvec ind = sort_index(real(eval));
    //std::vector< std::complex<double> > unique_eval={ eval(ind(0)) };
    vector<EigenStateC> res;
    for(uint i = 0; i < eval.size(); i++)
        res.push_back(EigenStateC {real(eval(ind(i))),evec.col(ind(i)),nPart,b.sym} );

    return  res;
}




inline EigenStateC FindGS(QOperatorC& ham, const FockBasis& b,int nPart, bool add = false)
{
    auto H = ham.toMatrix<std::complex<double> >(b);
    //H = 0.5*(H + H.t());
    cout<<"matrix nonzero="<<H.n_nonzero<<" Hermitian(tol=1e-4): "<<H.is_hermitian(1e-4)<<endl; cout.flush();

    cx_vec eval;
    cx_mat evec;

    eigs_gen(eval,evec, H, 1,"sr");

    EigenStateC res(real(eval(0)),evec.col(0),nPart,b.sym);
    return  res;
}



inline void PrintState(const Col<double> &state,const FockBasis &b,double tol=1e-3)
{

    uvec id=sort_index(abs(state),"descend");
    cout<<"State  Probability Coeff"<<endl;
    for(size_t i=0;i<state.size();i++)
        if(fabs(state(id(i)))>tol)
        {
            auto p=id(i);
            cout<<b.vec[p]<<" "<< pow(fabs(state(p)),2)<<" "<<state(p)<<endl;
        }
}



inline void PrintStateLLform(const Col<double> &state,const FockBasis &b, int Nl, double tol=1e-3)
{
    int L = b.vec.size();
    uvec id=sort_index(abs(state),"descend");

    for(size_t i=0;i<state.size();i++)
        if(fabs(state(id(i)))>tol)
        {
            auto p=id(i);
            cout<<"\nProbability, Coeff: ";
            cout<<pow(fabs(state(p)),2)<<", "<<state(p)<<endl;
            for(int j = L-2; j >= 0; j--)
            {

                cout<<b.vec[p][j];
                if(j%Nl==0) cout<<endl;

            }
            cout<<"-Eg-"<<endl;
            cout<<"000"<<b.vec[p][L-1]<<endl;
            cout<<endl<<endl;

        }
}


struct TimeEvolutionH
{
    mat evec;
    vec eval;
    vec psi0_n; //to calculate <n|psi0>

    TimeEvolutionH(const sp_mat &H, const FockBasis& b, const vec &psi0)
        :psi0_n(psi0.size())
    {
        //        auto H = ham.toMatrix<L>(b);
        cout<<"matrix nonzero="<<H.n_nonzero<<endl;
        mat evec1;
        vec eval1;

        eigs_sym(eval, evec, H , psi0.size()-1, "sa");
        eigs_sym(eval1, evec1, H , 1, "la");
        eval.insert_rows(psi0.size()-1, eval1);
        evec.insert_cols(psi0.size()-1,evec1.col(0));

        //eig_sym(eval,evec, mat(H));
        for(uint n=0;n<eval.size();n++)
            psi0_n[n]=dot(evec.col(n),psi0);

        //        for(uint n=0;n<eval.size();n++)
        //           cout<<eval1(n)<<" "<<eval(n)<<endl;


    }

    cx_vec EvolveTo(double t,double e0 = 0) const
    {
        cx_vec gt(eval.size());
        for(uint n=0;n<eval.size();n++)
        {
            cmpx q={0,(e0-eval[n])*t};
            gt(n)=psi0_n[n]*exp(q);
        }
        return evec*gt;
    }
};


static  cx_vec EvolveTridiagonal(cx_vec &a, cx_vec &b, double dt, double e_ref)
{
    int size = a.size();
    //cx_mat h(size,size,fill::zeros);
    sp_mat h(size,size);
    h.diag() = real(a);
    if (size>=2)
    {
        h.diag(-1) = real(b);
        h.diag(1) = real(b);
    }
    mat evec, evec1;
    vec eval, eval1;

    eigs_sym(eval, evec, h , size-1, "sa");
    eigs_sym(eval1, evec1, h, 1, "la");
    eval.insert_rows(size-1, eval1);
    evec.insert_cols(size-1,evec1.col(0));

    // eig_sym(eval, evec, h);

    //    cx_vec psi0(size);
    //    psi0(0) = 1;
    //    cx_vec psi0_n = evec.t()*psi0; //to calculate <n|psi0>
    vec psi0_n = evec.row(0).t(); //to calculate <n|psi0>


    cx_vec gt(eval.size());
    for(uint n=0; n<eval.size(); n++)
    {
        cmpx q = {0,(e_ref-eval[n])*dt};
        gt(n)=psi0_n[n]*exp(q);
    }
    return evec*gt;

    //   cmpx q = {0,-dt};
    //   if(!(h.is_hermitian(1e-20)))
    //   (h.t()-h).print("M");
    //   return expmat_sym(q*h)*psi0;
}

static cx_vec KrylovEvolveTo(const sp_mat &H, double dt, cx_vec &psi0, double e_ref = 0)
{
    int Nmax = std::min(15, int(psi0.size()-1));
    std::vector<cx_vec> v;
    cx_vec c;
    cx_vec a(Nmax), b(Nmax-1);

    double norm0 = Norm(psi0);
    v.push_back(psi0/norm0);
    a(0) = cdot(v[0],H*v[0]);
    for(auto j = 1; j < Nmax; j++)
    {
        cx_vec w = H*v[j-1];
        Orthogonalize(w,v,j);
        v.push_back(w/Norm(w));
        a(j) = cdot(v[j],H*v[j]);
        b(j-1)= cdot(v[j],H*v[j-1]);
        ///test convergence here
    }
    //    a.print("a");
    //    b.print("b");
    c = EvolveTridiagonal(a, b, dt, e_ref);
    cx_vec res(psi0.size());
    for(uint i = 0; i < c.size(); i++)
        res += c(i)*v[i];

    return norm0*res;
}

#endif // HAMILTONIANDIAGONALIZATION_H
