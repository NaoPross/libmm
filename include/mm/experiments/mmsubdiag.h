#pragma once

#ifndef __TRIDIAG_H__
#define __TRIDIAG_H__

template<class T>
class tridiag
{
    std::vector<T> diag, lower, upper;

    struct row
    {
        std::size_t index;
        tridiag& ref;

        T shared_zero = 0;

        row(std::size_t i, tridiag& _ref) : index(i), ref(_ref) {}

        T& operator[](std::size_t j)
        {
            switch(static_cast<long int>(j) - index)
            {
            case -1:
                return ref.lower[j];
            case 0:
                return ref.diag[j];
            case 1:
                return ref.upper[index];
            default:
                shared_zero = 0; // assure zero
                return shared_zero;
            }
        }
    };

    struct const_row
    {
        std::size_t index;
        const tridiag& ref;

        const T shared_zero = 0;

        const_row(std::size_t i, const tridiag& _ref) : index(i), ref(_ref) {}

        const T& operator[](std::size_t j) const
        {
            switch(static_cast<long int>(j) - index)
            {
            case -1:
                return ref.lower[j];
            case 0:
                return ref.diag[j];
            case 1:
                return ref.upper[index];
            default:
                return shared_zero;
            }
        }
    };

public:

    tridiag(std::size_t N = 1) : diag(N), lower(N-1), upper(N-1)
    {
    } 

    void resize(std::size_t N)
    {
        diag.reserve(N);
        lower.reserve(N-1);
        upper.reserve(N-1);
    }

    void zeros()
    {
        for (std::size_t i = 0; i < size()-1; ++i)
                diag[i] = lower[i] = upper[i] = 0;

        diag[size()-1] = 0;
    }

    // set diagonal to zero
    void zeros_diag()
    {
        for (std::size_t i = 0; i < size(); ++i)
            diag[i] = 0;
    }

    row operator[](std::size_t i)
    {
        return row(i, *this);
    }

    const const_row operator[](std::size_t i) const
    {
        return const_row(i, *this);
    }

    std::size_t size() const
    {
        return diag.size();
    }

    template<template <typename> class V>
    V<T> solve(const V<T>& rhs)
    {
        V<T> solution(diag.size());
        V<T> new_diag(diag);
        V<T> new_rhs(rhs);

        for(std::size_t i(1); i < size(); ++i)
        {
            T pivot = lower[i-1]/new_diag[i-1];
            new_diag[i] -= pivot * upper[i-1];
            new_rhs[i] -= pivot * new_rhs[i-1];
        }

        solution[size()-1] = new_rhs[size()-1] / new_diag[size()-1];

        for(int i(static_cast<int>(size()-2)); i>=0; --i)
            solution[i] = (new_rhs[i] - upper[i]*solution[i+1]) / new_diag[i];

        return solution;
    }

    static constexpr tridiag id(size_t N)
    {
        tridiag<T> diag;

        for (size_t k = 0; k < N; ++k)
            diag[k][k] = 1;

        return diag;
    }

    const std::vector<T>& d() const
    {
        return diag;
    }

    const std::vector<T>& down() const
    {
        return lower;
    }

    const std::vector<T>& up() const
    {
        return upper;
    }
};

// multiplication operator overloading
/*template<class T, template <class> class V>
V<T> operator*(const tridiag<T>& M, const V<T>& v)
{
    const std::size_t N = v.size();

    V<T> u(N);

    u[0] = M[0][0] * v[0] + M[0][1] * v[1];

    for (std::size_t k = 1; k < N-1; ++k)
        u[k] = M[k][k-1] * v[k-1] + M[k][k] * v[k] + M[k][k+1] * v[k+1];

    u[N-1] = M[N-1][N-2] * v[N-2] + M[N-1][N-1] * v[N-1];

    return u;
}*/

template<class T, template <class> class V>
V<T> operator*(const tridiag<T>& M, const V<T>& v)
{
    const std::size_t N = v.size();

    V<T> u(N);

    u[0] = M.d()[0] * v[0] + M.up()[0] * v[1];

    for (std::size_t k = 1; k < N-1; ++k)
        u[k] = M.down()[k] * v[k-1] + M.d()[k] * v[k] + M.up()[k] * v[k+1];

    u[N-1] = M.down()[N-1] * v[N-2] + M.d()[N-1] * v[N-1];

    return u;
}

template<class T>
tridiag<T> operator+(tridiag<T> M, const tridiag<T>& A)
{
    if (M.size() != A.size())
        return M;

    const std::size_t N = M.size();

    M[0][0] += A[0][0];
    M[0][1] += A[0][1];

    for (std::size_t k = 1; k < N-1; ++k)
        for (std::size_t j = k-1; j <= k+1; ++j)
            M[k][j] += A[k][j];

    M[N-1][N-2] += A[N-1][N-2];
    M[N-1][N-1] += A[N-1][N-1];

    return M;
}

template<class T>
tridiag<T> operator*(tridiag<T> M, const T& value)
{
    const std::size_t N = M.size();

    M[0][0] *= value;
    M[0][1] *= value;

    for (std::size_t k = 1; k < N-1; ++k)
        for (std::size_t j = k-1; j <= k+1; ++j)
            M[k][j] *= value;

    M[N-1][N-2] *= value;
    M[N-1][N-1] *= value;

    return M;
}

// add value * Id
template<class T>
tridiag<T> operator+(tridiag<T> M, const T& value)
{
    const std::size_t N = M.size();

    for (std::size_t k = 0; k < N; ++k)
        M[k][k] += value;

    return M;
}

template<class T>
tridiag<T> operator*(const T& value, const tridiag<T>& M)
{
    return M * value;
}

// add value * Id
template<class T>
tridiag<T> operator+(const T& value, const tridiag<T>& M)
{    
    return M + value;
}

#include <ostream>

template<class T>
std::ostream& operator<<(std::ostream& os, const tridiag<T>& M)
{
    const std::size_t N = M.size();


    os << M[0][0] << " ";
    os << M[0][1] << " ";
    os << std::endl;

    for (std::size_t k = 1; k < N-1; ++k) {
        for (std::size_t j = k-1; j <= k+1; ++j)
            os << M[k][j] << " ";
        os << std::endl;
    }


    os << M[N-1][N-2] << " ";
    os << M[N-1][N-1] << std::endl;

    return os;
}

#endif
