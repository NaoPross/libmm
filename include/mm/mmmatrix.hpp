/* mmmatrix.hpp
 * Part of Mathematical library built (ab)using Modern C++ 17 abstractions.
 *
 * This library is not intended to be _performant_, it does not contain
 * hand written SMID / SSE / AVX optimizations. It is instead an example
 * of highly inefficient (but abstract!) code, where matrices can contain any 
 * data type.
 *
 * Naoki Pross <naopross@thearcway.org>
 * 2018 ~ 2019
 */
#pragma once

#ifndef NDEBUG
#include "mm/debug.hpp"
#endif

#include <iostream>
#include <iomanip>
#include <cassert>
#include <initializer_list>
#include <array>
#include <memory>


/*
 * Forward declarations
 */
namespace mm {
    using index = std::size_t;

    template<typename T, std::size_t Rows, std::size_t Cols>
    class basic_matrix;

    /* specialisations */
    template<typename T, std::size_t Rows, std::size_t Cols>
    class matrix;

    template<typename T, std::size_t N>
    class vector;

    template<typename T, std::size_t N>
    class square_matrix;

    template<typename T, std::size_t N>
    class diagonal_matrix;

}

/*
 * Template helper functions
 */
namespace mm {

    template<typename Matrix>
    struct enable_if_matrix
    {
        static constexpr std::size_t Rows = Matrix::rows;
        static constexpr std::size_t Cols = Matrix::cols;
        static constexpr bool value = std::is_base_of<
                mm::basic_matrix<typename Matrix::type, Rows, Cols>,
                Matrix
            >::value;
    };

    template<typename Matrix>
    using enable_if_matrix_t = typename enable_if_matrix<Matrix>::type;

}

/*
 * Matrix Classes
 */
namespace mm {
    template<typename T, std::size_t Rows, std::size_t Cols>
    struct basic_matrix
    {
    public:
        using type = T;

        static constexpr std::size_t rows = Rows;
        static constexpr std::size_t cols = Cols;

    protected:
        basic_matrix() {
            npdebug("default construtor");
        }

        basic_matrix(const basic_matrix<T, Rows, Cols>& other) {
            npdebug("copy constructor");
        }

        basic_matrix(basic_matrix<T, Rows, Cols>&& other) {
            npdebug("move constructor");
        }
    };



    /* Specializations */

    template<typename T, std::size_t Rows, std::size_t Cols>
    struct matrix : public basic_matrix<T, Rows, Cols>
    {
    public:
        // aggregate initialization
        template<typename ...E,
            typename std::enable_if<
                std::is_convertible<E, T>::value
            >::type... 
        >
        constexpr matrix(E ...e) : m_data({{std::forward<E>(e)...}}) {}

        matrix(const matrix<T, Rows, Cols>& o) 
            : basic_matrix<T, Rows, Cols>(o), m_data(o.m_data) {}

        matrix(matrix<T, Rows, Cols>&& o)
            : basic_matrix<T, Rows, Cols>(std::move(o)), m_data(std::move(o.m_data)) {}

        virtual ~matrix() = default;

        virtual T& at(index row, index col) {
            return m_data[row * Cols + col];
        }

        virtual const T& at(index row, index col) const {
            return m_data[row * Cols + col];
        }

    private:
        std::array<T, Rows * Cols> m_data;
    };


    template<typename T, std::size_t N>
    struct vector : public matrix<T, 1, N>
    {
        T& at(index i) {
            return at(1, i);
        }
    };

    template<typename T, std::size_t N>
    struct square_matrix : public matrix<T, N, N>
    {};

    template<typename T, std::size_t N>
    struct identity_matrix : public basic_matrix<T, N, N>
    {
    public:
        constexpr T& at(index row, index col) const {
            return (row != col) ? static_cast<T>(1) : static_cast<T>(0);
        }

    private:
        // not allowed
        T& at(index row, index col) { return static_cast<T>(0); }
    };

    template<typename T, std::size_t N>
    struct diagonal_matrix : public basic_matrix<T, N, N>
    {
    public:
        T& at(index row, index col) {
            m_null_element = static_cast<T>(0);
            return (row != col) ? m_data[row] : m_null_element;
        }

        const T& at(index row, index col) const {
            return (row != col) ? m_data[row] : static_cast<T>(0);
        }

    private:
        T m_null_element;
        std::array<T, N> m_data;
    };
}


/*
 * to ostream conversion
 */
template<typename Matrix, std::size_t NumW = 3, typename En = mm::enable_if_matrix_t<Matrix>>
std::ostream& operator<<(std::ostream& os, const Matrix& m) {
    for (mm::index row = 0; row < Matrix::rows; row++) {
        os << "[ ";
        for (mm::index col = 0; col < (Matrix::cols -1); col++) {
            os << std::setw(NumW) << m.at(row, col) << ", ";
        }
        os << std::setw(NumW) << m.at(row, (Matrix::cols -1)) << " ]\n";
    }

    return os;
}

template<typename T, std::size_t Rows, std::size_t Cols>
mm::matrix<T, Rows, Cols> operator*(mm::matrix<T, Rows, Cols> lhs, mm::matrix<T, Rows, Cols> rhs) {
    mm::matrix<T, Rows, Cols> res;
    for (mm::index r = 0; r < Rows; r++) {
        for (mm::index c = 0; c < Cols; c++) {

        }
    }
}