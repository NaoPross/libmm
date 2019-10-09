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

#include "mm/debug.hpp"

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


        template<typename U, std::size_t ORows, std::size_t OCols>
        friend class mm::basic_matrix;

        virtual ~basic_matrix() {};

        // copy from another matrix
        // template<std::std::size_t ORows, std::size_t OCols>
        // basic_matrix(const basic_matrix<T, ORows, OCols>& other) {
        //     static_assert(ORows <= Rows);
        //     static_assert(OCols <= Cols);

        //     for (index row = 0; row < Rows; row++)
        //         for (index col = 0; col < Cols; col++)
        //             at(row, col) = other.at(row, col);
        // }

        virtual T& at(index row, index col) = 0;
        virtual const T& at(index row, index col) const  = 0;

        // constexpr std::size_t rows() { return Rows; }
        // constexpr std::size_t cols() { return Cols; }

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
        matrix(E ...e) : m_data({{std::forward<E>(e)...}}) {}

        matrix(const matrix<T, Rows, Cols>& o) 
            : basic_matrix<T, Rows, Cols>(o), m_data(o.m_data) {}

        matrix(matrix<T, Rows, Cols>&& o)
            : basic_matrix<T, Rows, Cols>(std::move(o)), m_data(std::move(o.m_data)) {}

        virtual ~matrix() = default;

        virtual T& at(index row, index col) override {
            return m_data[row * Cols + col];
        }

        virtual const T& at(index row, index col) const override {
            return m_data[row * Cols + col];
        }

    private:
        std::array<T, Rows * Cols> m_data;
    };


    template<typename T, std::size_t N>
    struct vector : public matrix<T, 1, N> {};

    template<typename T, std::size_t N>
    struct square_matrix : public matrix<T, N, N>
    {};

    template<typename T, std::size_t N>
    struct identity_matrix : public basic_matrix<T, N, N>
    {
    public:
        const T& at(index row, index col) const override {
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
        T& at(index row, index col) override {
            m_null_element = static_cast<T>(0);
            return (row != col) ? m_data[row] : m_null_element;
        }

        const T& at(index row, index col) const override {
            return (row != col) ? m_data[row] : static_cast<T>(0);
        }

    private:
        T m_null_element;
        std::array<T, N> m_data;
    };
}

/*
 *  Matrix Opertors
 */

namespace mm {
}


template<typename T, std::size_t Rows, std::size_t Cols, unsigned NumW = 3>
std::ostream& operator<<(std::ostream& os, const mm::basic_matrix<T, Rows, Cols>& m) {
    for (mm::index row = 0; row < Rows; row++) {
        os << "[ ";
        for (mm::index col = 0; col < (Cols -1); col++) {
            os << std::setw(NumW) << m.at(row, col) << ", ";
        }
        os << std::setw(NumW) << m.at(row, (Cols -1)) << " ]\n";
    }

    return os;
}
   