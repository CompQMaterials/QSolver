#ifndef HAMSOLVER_H
#define HAMSOLVER_H

#include<iostream>
#include<iomanip>
#include<vector>
#include<complex>
#include<fstream>
#include<string>
#include<armadillo>

#include"fockbasis.h"
#include"qoperator.h"
#include"symmetrygroup.h"
#include"lanczos.h"


template<class scalar=double>
struct EigenStateG
{
    double ener=std::numeric_limits<double>::max();
    Col<scalar> state;
    int nPart;
    int sym;

    EigenStateG(){}
    EigenStateG(double ener,const Col<scalar>& state,int nPart,int sym)
        :ener(ener),state(state), nPart(nPart), sym(sym) {}

    void Save(std::string filename) const
    {
        std::ofstream out(filename);
        out<<std::setprecision(14);
        out<<ener<<" ";
        out<<nPart<<" ";
        out<<sym<<endl;
        state.save(out);
    }
    bool Load(std::string filename)
    {
        std::ifstream in(filename);
        if (!in) return false;
        in>>ener;
        in>>nPart;
        in>>sym;
        state.load(in, arma::arma_binary);
        return true;
    }

    void PrintToFile(const FockBasis &b, std::ofstream &out,double tol=1e-3) const
    {

        uvec id=sort_index(abs(state),"descend");

        out<<"Energy:"<<ener<<endl;
        out<<"State  Probability Coeff"<<endl;
        for(size_t i=0;i<state.size();i++)
            if(fabs(state(id(i)))>tol)
            {
                auto p=id(i);
                out<<b.vec[p]<<" "<< pow(fabs(state(p)),2)<<" "<<state(p)<<endl;
            }
    }

    void Print(const FockBasis &b,double tol=1e-3) const
    {

        uvec id=sort_index(abs(state),"descend");
        cout<<std::setprecision(10);
        cout<<"Energy:"<<ener<<endl;
        cout<<"State  Probability Coeff"<<endl;
        for(size_t i=0;i<state.size();i++)
            if(fabs(state(id(i)))>tol)
            {
                auto p=id(i);
                cout<<b.vec[p]<<" "<< pow(fabs(state(p)),2)<<" "<<state(p)<<endl;
            }
    }
};

template<class scalar>
bool operator<(const EigenStateG<scalar>& e1,const EigenStateG<scalar>& e2)
{
    return e1.ener<e2.ener;
}

using EigenState=EigenStateG<>;


inline EigenState FindGS(QOperator ham,FockBasisFixedChargeG& b,
                  const SymmetryGroup<double>& G, int nPart)
{
    auto H=ham.template toMatrix(b,G);
    cout<<"matrix nonzero="<<H.n_nonzero<<endl; cout.flush();

    vec eval;
    mat evec;
    if (b.Size()>10)
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
        if (unique_eval.size()>10) break; // only print the first 10 different evals
        cout<<x<<"\n";
    }
    cout<<"dim="<<b.Size()<<" nPart="<<nPart<<" sym="<<b.sym<<" ener="<<eval(ind(0))<<endl;
    return {eval(ind(0)),evec.col(ind(0)),nPart,b.sym};

    //    Col<double> wf(b.Size());
    //    wf.randu();
    //    auto sol=Diagonalize<double>(H,wf);
    //    cout<<sol.cIter<<" lanczos iter; "<<"nPart="<<nPart<<" sym="<<b.sym<<" ener="<<sol.ener0<<endl;
    //    return {sol.ener0,sol.x0,nPart,b.sym};
}


inline EigenStateG<cmpx> FindGS(QOperatorG<cmpx> ham,FockBasisFixedChargeG& b,
                         const SymmetryGroup<cmpx>& G, int nPart)
{
    //        b.SetSym(G,sym);
    auto H=ham.toMatrix<cmpx>(b,G);
    cout<<"matrix nonzero="<<H.n_nonzero<<endl; cout.flush();

    cx_vec evalz;
    vec eval;
    cx_mat evec;
    if(b.Size()>10)
    {
        eigs_gen(evalz,evec,H ,std::min(b.Size()-2,201),"sr");
        eval=arma::real(evalz);
    }
    else
        eig_sym(eval,evec,cx_mat(H));
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
        if (unique_eval.size()>10) break; // only print the first 10 different evals
        cout<<x<<"\n";
    }
    cout<<"dim="<<b.Size()<<" nPart="<<nPart<<" sym="<<b.sym<<" ener="<<eval(ind(0))<<endl;
    return {eval(ind(0)),evec.col(ind(0)),nPart,b.sym};

    //    Col<cmpx> wf(b.Size());
    //    wf.randu();
    //    auto sol=Diagonalize<cmpx>(H,wf);
    //    cout<<sol.cIter<<" lanczos iter; "<<"nPart="<<nPart<<" sym="<<b.sym<<" ener="<<sol.ener0<<endl;
    //    return {sol.ener0,sol.x0,nPart,b.sym};
}



inline EigenStateG<cmpx> FindGS(QOperator ham,FockBasisFixedChargeG& b,
                         const SymmetryGroup<cmpx>& G, int nPart)
{
    //        b.SetSym(G,sym);
    auto H=ham.toMatrix<cmpx>(b,G);
    cout<<"matrix nonzero="<<H.n_nonzero<<endl; cout.flush();

    cx_vec evalz;
    vec eval;
    cx_mat evec;
    if(b.Size()>10)
    {
        eigs_gen(evalz,evec,H ,std::min(b.Size()-2,101),"sr");
        eval=arma::real(evalz);
    }
    else
        eig_sym(eval,evec,cx_mat(H));
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
    return {eval(ind(0)),evec.col(ind(0)),nPart,b.sym};

    //    Col<cmpx> wf(b.Size());
    //    wf.randu();
    //    auto sol=Diagonalize<cmpx>(H,wf);
    //    cout<<sol.cIter<<" lanczos iter; "<<"nPart="<<nPart<<" sym="<<b.sym<<" ener="<<sol.ener0<<endl;
    //    return {sol.ener0,sol.x0,nPart,b.sym};
}



inline EigenState FindGS(QOperator ham,FockBasisFixedCharge& b,int nPart, int nstates = 101)
{
    auto H=ham.toMatrix(b);
    cout<<"matrix nonzero="<<H.n_nonzero<<endl; cout.flush();

    vec eval;
    mat evec;
    if(b.Size()>10)
        eigs_sym(eval,evec, H,std::min(b.Size()-1,nstates),"sa");
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
//        if (unique_eval.size()>10) break; // only print the first 2 different evals
//        cout<<x/*-eval(ind(0))*/<<"\n";
    }
    cout<<"dim="<<b.Size()<<" nPart="<<nPart<<" sym="<<b.sym<<" ener="<<eval(ind(0))<<endl;
    return {eval(ind(0)),evec.col(ind(0)),nPart,b.sym};

    //    Col<double> wf(b.Size());
    //    wf.randu();
    //    auto sol=Diagonalize<double>(H,wf);
    //    cout<<sol.cIter<<" lanczos iter; "<<"nPart="<<nPart<<" sym="<<b.sym<<" ener="<<sol.ener0<<endl;
    //    return {sol.ener0,sol.x0,nPart,b.sym};
}




//----------------------------------------- Time evolution---------------------


struct TimeEvolution
{
    cx_mat evec;
    vec eval;
    cx_vec psi0_n; //to calculate <n|psi0>

    TimeEvolution(const sp_cx_mat &H, const cx_vec &psi0)
        :psi0_n(psi0.size())
    {
        eig_sym(eval,evec,cx_mat(H));
        for(uint n=0;n<eval.size();n++)
            psi0_n[n]=cdot(evec.col(n),psi0);
    }
    cx_vec EvolveTo(double t) const
    {
        cx_vec gt(eval.size());
        for(uint n=0;n<eval.size();n++)
        {
            cmpx q={0,-eval[n]*t};
            gt(n)=psi0_n[n]*exp(q);
        }
        return evec*gt;
    }
};



#endif // HAMSOLVER_H
