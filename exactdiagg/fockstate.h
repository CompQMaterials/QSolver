#ifndef FOCKSTATE_H
#define FOCKSTATE_H
#include<vector>
#include<map>
#include<unordered_map>
#include<iostream>
#include<limits>
#include<functional>
#include<cmath>
#include<complex>
#include <boost/dynamic_bitset.hpp>

using FockState = boost::dynamic_bitset<>;

inline FockState TensorProd(const FockState& f1, const FockState& f2)
{
    FockState f1p = f1;
    FockState f2p = f2;
    f1p.resize(f1.size() + f2.size()); // Ensure f1p has enough space
    return (f1p << f2.size()) | f2p;
}

inline std::vector<FockState> Split(const FockState& f, int R)
{
    const int step = f.size() / R;
    FockState mask((boost::dynamic_bitset<>::size_type(1) << step) - 1);

    std::vector<FockState> a(R);
    FockState temp = f;
    for(int i = 0; i < R; i++)
    {
        a[i] = temp & mask;
        temp >>= step;
    }
    return a;
}

inline FockState Join(const std::vector<FockState>& a)
{
    int L1 = a[0].size();
    FockState f(L1 * a.size());

    for(int i = 0; i < a.size(); i++)
    {
        FockState mask = a[i];
        mask.resize(L1 * a.size()); // Ensure mask has enough space
        f |= (mask << (i * L1));
    }
    return f;
}

inline bool BitParityLeft(FockState f, int j)
{
    f >>= (j + 1);
    return f.count() & 1;
}

struct FockLess {
    bool operator()(const FockState& f1, const FockState& f2) const
    {
        for(int i = f1.size() - 1; i > 0; i--)
        {
            if(f1[i] != f2[i]) return f1[i] < f2[i];
        }
        return false;
    }
};

inline bool operator<(const FockState& f1, const FockState& f2)
{
    return f1.to_ulong() < f2.to_ulong();
}

template<class scalar>
using State = std::map<FockState, scalar, FockLess>;

template<class scalar>
double Norm(const State<scalar>& s)
{
    double sum = 0;
    for(const auto& it : s)
        sum += std::norm(it.second);
    return std::sqrt(sum);
}

#endif // FOCKSTATE_H
