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

#include <iostream>
#include <cassert>
#include <initializer_list>
#include <array>
#include <memory>


namespace mm {
    using index = std::size_t;

    template<typename T, std::size_t Rows, std::size_t Cols>
    class basic_matrix;

    /* specialisations */
    template<typename T, std::size_t Rows, std::size_t Cols>
    class matrix;

    template<typename T, std::size_t N>
    class square_matrix;

    template<typename T, std::size_t N>
    class diagonal_matrix;
}

/*
 * Matrix class, no access methods
 */
namespace mm {
    template<typename T, std::size_t Rows, std::size_t Cols>
    class basic_matrix
    {
    public:
        using type = T;

        template<typename U, std::size_t ORows, std::size_t OCols>
        friend class mm::matrix;

        // copy from another matrix
        template<std::size_t ORows, std::size_t OCols>
        matrix(const basic_matrix<T, ORows, OCols>& other);

        virtual T& at(index row, index col) = 0;
        virtual const T& at(index row, index col) const = 0;
    };


    /* Specializations */

    template<typename T, std::size_t Rows, std::size_t Cols>
    struct matrix : public basic_matrix<T, N>
    {
    public:
        virtual T& at(index row, index col) override {
            return m_data[row * Cols + col];
        }

        virtual const T& at(index row, index col) const override {
            return at(row, col);
        }

    private:
        std::array<T, Rows * Cols> m_data;
    };


    template<typename T, std::size_t N>
    struct vector : public matrix<T, 1, N> {};

    template<typename T, std::size_t N>
    struct square_matrix : public basic_matrix<T, N>
    {
    public:
        virtual T& at(index row, index col) override {
            return m_data[row * N + col];
        }

        virtual const T& at(index row, index col) const override {
            return at(row, col);
        }

    private:
        std::array<T, N*N> m_data;
    };

    template<typename T, std::size_t N>
    struct identity_matrix : public basic_matrix<T, N, N>
    {
    public:
        const T& at(index row, index col) const override {
            return (row != col) ? static_cast<T>(1) : static_cast<T>(0);
        }

    private:
        T m_useless;
        T& at(index row, index col) { return m_useless; }
    }

    template<typename T, std::size_t N>
    struct diagonal_matrix : public basic_matrix<T, N, N>
    {
    public:
        T& at(index row, index col) override {
            n_null_element = static_cast<T>(0);
            return (row != col) ? m_data[row] : n_null_element;
        }

        const T& at(index row, index col) const override {
            return (row != col) ? m_data[row] : static_cast<T>(0);
        }

    private:
        T m_null_element;
        std::array<T, N> m_data;
    }
}
