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
#include <iomanip>
#include <cstring>
#include <cassert>
#include <initializer_list>
#include <array>
#include <memory>

#include "mm/mmiterator.hpp"

namespace mm {

    template<typename T, std::size_t Rows, std::size_t Cols>
    class basic_matrix;

    /* specialisations */

    template<typename T, std::size_t Rows, std::size_t Cols>
    class matrix; // simple matrix format

    /* specialization of basic_matrx for Cols = 1 */
    template<typename T, std::size_t Rows>
    class row_vec;

    /* specialization of basic_matrx for Rows = 1 */
    template<typename T, std::size_t Cols>
    class col_vec; // transposed version of row_vec

    template<typename T, std::size_t N>
    class square_matrix;

    /* specialisation of a square_matrix for a sub-diagonal composed matrix */
    template<typename T, std::size_t N, std::size_t K = 0>
    class diag_matrix;
}


/*namespace mm {

    template<typename T, std::size_t N>
    using diag_iterator = vector_iterator<T, N, N, MM_DIAG_ITER, mm::basic_matrix<T, N, N>>;

    template<typename T, std::size_t N>
    using const_diag_iterator = vector_iterator<typename std::add_const<T>::type, N, N, MM_DIAG_ITER, typename std::add_const<mm::basic_matrix<T, N, N>>::type>;
}*/


/*
 * Matrix class, no access methods
 */

template<typename T, std::size_t Rows, std::size_t Cols>
class mm::basic_matrix
{
public:
    using type = T;

    template<typename U, std::size_t ORows, std::size_t OCols>
    friend class mm::basic_matrix;

    template<typename U, std::size_t ORows, std::size_t OCols>
    friend class mm::matrix;

    template<typename U, std::size_t ORows, std::size_t OCols, class Grid>
    friend class mm::iter::basic_iterator;

    template<typename U, std::size_t ON, class Grid>
    friend class mm::iter::diag_iterator;

    //template<typename U, std::size_t ORows, std::size_t OCols, class Grid>
    //friend class mm::iter::basic_iterator<T, Rows, Cols, mm::basic_matrix<T, Rows, Cols>>;

    //template<typename U, std::size_t ORows, std::size_t OCols, class Grid>
    //friend class mm::iter::basic_iterator<typename std::add_const<T>::type, Rows, Cols, typename std::add_const<mm::basic_matrix<T, Rows, Cols>>::type>;

    basic_matrix();

    // from initializer_list
    basic_matrix(std::initializer_list<std::initializer_list<T>> l);

    // copyable and movable
    basic_matrix(const basic_matrix<T, Rows, Cols>& other) = default;
    basic_matrix(basic_matrix<T, Rows, Cols>&& other) = default;

    // copy from another matrix
    template<std::size_t ORows, std::size_t OCols>
    basic_matrix(const basic_matrix<T, ORows, OCols>& other);

    void swap_rows(std::size_t x, std::size_t y);
    void swap_cols(std::size_t x, std::size_t y);

    // mathematical operations
    //virtual basic_matrix<T, Cols, Rows> transposed() const;
    //inline basic_matrix<T, Cols, Rows> td() const { return transposed(); }

protected:
    template<typename ConstIterator>
    basic_matrix(ConstIterator begin, ConstIterator end);

private:
    std::array<T, Rows * Cols> data;
};


template<typename T, std::size_t Rows, std::size_t Cols>
mm::basic_matrix<T, Rows, Cols>::basic_matrix() {
    std::fill(data.begin(), data.end(), 0);
}

template<typename T, std::size_t Rows, std::size_t Cols>
mm::basic_matrix<T, Rows, Cols>::basic_matrix(
    std::initializer_list<std::initializer_list<T>> l
) {
    assert(l.size() == Rows);
    auto data_it = data.begin();

    for (auto&& row : l) {
        data_it = std::copy(row.begin(), row.end(), data_it);
    }
}

template<typename T, std::size_t Rows, std::size_t Cols>
template<std::size_t ORows, std::size_t OCols>
mm::basic_matrix<T, Rows, Cols>::basic_matrix(
    const mm::basic_matrix<T, ORows, OCols>& other
) {
    static_assert((ORows <= Rows),
        "cannot copy a taller matrix into a smaller one"
    );

    static_assert((OCols <= Cols),
        "cannot copy a larger matrix into a smaller one"
    );

    std::fill(data.begin(), data.end(), 0);
    for (unsigned row = 0; row < Rows; row++)
        for (unsigned col = 0; col < Cols; col++)
            this->at(row, col) = other.at(row, col);
}

/* protected construtor */
template<typename T, std::size_t Rows, std::size_t Cols>
template<typename ConstIterator>
mm::basic_matrix<T, Rows, Cols>::basic_matrix(
    ConstIterator begin, ConstIterator end
) {
    assert(static_cast<unsigned>(std::distance(begin, end)) >= ((Rows * Cols)));
    std::copy(begin, end, data.begin());
}

template<typename T, std::size_t Rows, std::size_t Cols>
void mm::basic_matrix<T, Rows, Cols>::swap_rows(std::size_t x, std::size_t y) {
    if (x == y)
        return;

    for (unsigned col = 0; col < Cols; col++)
        std::swap(this->at(x, col), this->at(y, col));
}

template<typename T, std::size_t Rows, std::size_t Cols>
void mm::basic_matrix<T, Rows, Cols>::swap_cols(std::size_t x, std::size_t y) {
    if (x == y)
        return;

    for (unsigned row = 0; row < Rows; row++)
        std::swap(this->at(row, x), this->at(row, y));
}

/*
 * Matrix object
 */

template<typename T, std::size_t Rows, std::size_t Cols>
class mm::matrix
{
public:

    //template<typename U, std::size_t ORows, std::size_t OCols>
    using vec_iterator = mm::iter::basic_iterator<T, Rows, Cols, mm::basic_matrix<T, Rows, Cols>>;

    //template<typename U, std::size_t ORows, std::size_t OCols>
    using const_vec_iterator = mm::iter::basic_iterator<typename std::add_const<T>::type, Rows, Cols, typename std::add_const<mm::basic_matrix<T, Rows, Cols>>::type>;
    
    // default zeros constructor
    matrix() : M(std::make_shared<mm::basic_matrix<T, Rows, Cols>>()), transposed(false) {}
    
    // from initializer_list
    matrix(std::initializer_list<std::initializer_list<T>> l)
        : M(std::make_shared<mm::basic_matrix<T, Rows, Cols>>(l)), transposed(false) {}

    // copyable and movable
    matrix(const matrix<T, Rows, Cols>& other) // deep copy
        : M(std::make_shared<mm::basic_matrix<T, Rows, Cols>>(*other.M)), transposed(other.transposed) {}

    matrix(basic_matrix<T, Rows, Cols>&& other) // move ptr
        : M(other.M), transposed(other.transposed)
    {
        other.M = nullptr;
    }

    // copy from another matrix
    /*template<std::size_t ORows, std::size_t OCols>
    matrix(const matrix<T, ORows, OCols>& other)
        : M(std::make_shared<mm::basic_matrix<T, Rows, Cols>(*other.M)), transposed(other.transposed) {} */

    matrix<T, Rows, Cols> operator=(const basic_matrix<T, Rows, Cols>& other) // deep copy
    {         
        *M = *other.M;
        transposed = other.transposed;
    }

    /* 
     * Transposition
     */

    matrix<T, Rows, Cols>& transpose_d()
    {
        transposed = !transposed;
        return *this;
    }

    matrix<T, Rows, Cols> transpose() const
    {
        auto m = shallow_cpy();
        m.transposed = !transposed;

        return m;
    }

    inline matrix<T, Rows, Cols>& td()
    {
        return transpose();
    }

    inline matrix<T, Rows, Cols> t() const
    {
        return transpose();
    }

    // strongly transpose
    matrix<T, Cols, Rows> transpose_cpy() const
    {
        matrix<T, Cols, Rows> out(); // copy
        // TODO
    }

    /*
     * Pointer status 
     */

    bool expired() const
    {
        return M == nullptr;
    }

    /*
     * Downcasting conditions
     */

    /// downcast to square matrix
    static inline constexpr bool is_square() { return (Rows == Cols); }
    inline constexpr square_matrix<T, Rows> to_square() const {
        static_assert(is_square());
        return static_cast<square_matrix<T, Rows>>(*this);
    }

    /// downcast to row_vector
    static inline constexpr bool is_row_vec() { return (Cols == 1); }
    inline constexpr row_vec<T, Rows> to_row_vec() const {
        static_assert(is_row_vec());
        return static_cast<row_vec<T, Rows>>(*this);
    }

    /// downcast to col_vector
    static inline constexpr bool is_col_vec() { return (Rows == 1); }
    inline constexpr col_vec<T, Cols> to_col_vec() const {
        static_assert(is_col_vec());
        return static_cast<col_vec<T, Cols>>(*this);
    }

    /* Accessors */

    T& at(std::size_t row, std::size_t col)
    {
        return (transposed) ? M->data[col * Cols + row] : M->data[row * Cols + col];
    }

    const T& at(std::size_t row, std::size_t col) const
    {
        return (transposed) ? M->data[col * Cols + row] : M->data[row * Cols + col];
    }

    std::size_t rows() const {
        return (transposed) ? Cols : Rows;
    }

    std::size_t cols() const {
        return (transposed) ? Rows : Cols;
    }

    mm::matrix<T, Rows, Cols>::vec_iterator operator[](std::size_t index)
    {
        return mm::matrix<T, Rows, Cols>::vec_iterator(*M, index, 0, !transposed);
    }

    mm::matrix<T, Rows, Cols>::const_vec_iterator operator[](std::size_t index) const
    {
        return mm::matrix<T, Rows, Cols>::const_vec_iterator(*M, index, 0, !transposed);
    }

    /*
     * Basic matematical operations (dimension indipendent)
     */

    mm::matrix<T, Rows, Cols>& operator+=(const mm::matrix<T, Rows, Cols>& m) {

        for (unsigned row = 0; row < std::min(rows(), m.rows()); ++row)
            for (unsigned col = 0; col < std::min(cols(), m.cols()); ++col)
                at(row, col) += m.at(row, col);
    
        return *this;
    }

    mm::matrix<T, Rows, Cols>& operator-=(const mm::matrix<T, Rows, Cols>& m) {

        for (unsigned row = 0; row < std::min(rows(), m.rows()); ++row)
            for (unsigned col = 0; col < std::min(cols(), m.cols()); ++col)
                at(row, col) -= m.at(row, col);
    
        return *this;
    }

    mm::matrix<T, Rows, Cols> operator*=(const T& k) {

        for (unsigned row = 0; row < rows(); ++row)
            for (auto& x : (*this)[row])
                x *= k;
    
        return *this;
    }

protected:

    std::shared_ptr<mm::basic_matrix<T, Rows, Cols>> M;

    // shallow construction
    matrix(std::shared_ptr<mm::basic_matrix<T, Rows, Cols>> grid, bool tr = false) : M(grid), transposed(tr) {}

    matrix<T, Rows, Cols> shallow_cpy() const
    {
        return matrix<T, Rows, Cols>(M, transposed);
    }

private:

    bool transposed;     
};

/* Basic operator overloading (dimension indipendent) */

template<typename T, std::size_t Rows, std::size_t Cols>
mm::matrix<T, Rows, Cols> operator+(
    mm::matrix<T, Rows, Cols> a,
    const mm::matrix<T, Rows, Cols>& b
) {
    return a += b;
}

template<typename T, std::size_t Rows, std::size_t Cols>
mm::matrix<T, Rows, Cols> operator-(
    mm::matrix<T, Rows, Cols> a,
    const mm::matrix<T, Rows, Cols>& b
) {
    return a -= b;
}

template<typename T, std::size_t Rows, std::size_t Cols>
mm::matrix<T, Rows, Cols> operator*(
    mm::matrix<T, Rows, Cols> a,
    const T& k
) {
    return a *= k;
}

template<typename T, std::size_t Rows, std::size_t Cols>
mm::matrix<T, Rows, Cols> operator*(
    const T& k,
    mm::matrix<T, Rows, Cols> a
) {
    return a *= k;
}

// simple multiplication
template<typename T, std::size_t M, std::size_t P1, std::size_t P2, std::size_t N>
mm::matrix<T, M, N> operator*(
    const mm::matrix<T, M, P1>& a,
    const mm::matrix<T, P2, N>& b
) {
    // TODO, adjust asserts for transposed cases
    static_assert(P1 == P2, "invalid matrix multiplication");
    assert(a.cols() == b.rows());

    mm::matrix<T, M, N> result;
    const mm::matrix<T, P2, N> bt = b.t(); // weak transposition

    npdebug("Calling *")

    for (unsigned row = 0; row < M; row++)
        for (unsigned col = 0; col < N; col++)
            result.at(row, col) = a[row] * bt[col]; // scalar product

    return result;
}

/*
 * Matrix operator <<
 */

template<typename T, std::size_t Rows, std::size_t Cols, unsigned NumW = 3>
std::ostream& operator<<(std::ostream& os, const mm::matrix<T, Rows, Cols>& m) {

    for (unsigned index = 0; index < m.rows(); index++) {
        os << "[ ";
        for (unsigned col = 0; col < m.cols()-1; ++col) {
            os << std::setw(NumW) << m.at(index, col) << ", ";
        }
        os << std::setw(NumW) << m.at(index, m.cols()-1) << " ]\n";
    }

    return os;
}

/*
 * Square matrix
 */

template<typename T, std::size_t N>
class mm::square_matrix : public mm::matrix<T, N, N>
{
public:

    using mm::matrix<T, N, N>::matrix;

    using diag_iterator = mm::iter::diag_iterator<T, N, mm::basic_matrix<T, N, N>>;

    using const_diag_iterator = mm::iter::diag_iterator<typename std::add_const<T>::type, N, typename std::add_const<mm::basic_matrix<T, N, N>>::type>;

    T trace();
    inline T tr() { return trace(); }

    mm::square_matrix<T, N>::diag_iterator diag_beg(int row = 0)
    {
        return diag_iterator(*(this->M), row, 0);
    }

    mm::square_matrix<T, N>::const_diag_iterator diag_end(int row = 0) const
    {
        return const_diag_iterator(*(this->M), row, N);
    }

    // TODO, determinant

    /// in place inverse
    // TODO, det != 0
    // TODO, use gauss jordan for invertible ones
    //void invert();, TODO, section algorithm
    
    /*
     * Generate the identity
     */

    static inline constexpr mm::square_matrix<T, N> identity() {
        mm::square_matrix<T, N> i;
        for (unsigned row = 0; row < N; row++)
            for (unsigned col = 0; col < N; col++)
                i.at(row, col) = (row == col) ? 1 : 0;

        return i;
    }
};

template<typename T, std::size_t N>
T mm::square_matrix<T, N>::trace() 
{
    T sum = 0;
    for (const auto& x : diag_beg())
        sum += x;

    return sum;
}


