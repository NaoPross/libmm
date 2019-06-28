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

namespace mm {

    template<typename T, std::size_t Rows, std::size_t Cols>
    class basic_matrix;

    template<typename T, std::size_t Rows, std::size_t Cols>
    class transposed_matrix;

    /* specialization of basic_matrx for Cols = 1 */
    template<typename T, std::size_t Rows>
    class row_vec;

    /* specialization of basic_matrx for Rows = 1 */
    template<typename T, std::size_t Cols>
    class col_vec;

    /* shorter name for basic_matrix */
    template<typename T, std::size_t Rows, std::size_t Cols>
    class matrix;

    /* specialization of basic_matrix for Rows == Cols */
    template<typename T, std::size_t N>
    class square_matrix;

    template<typename T, std::size_t N, int K = 0>
    class diagonal_matrix;

    /* 
     * Iterators 
     */

    /*
     * Most abstract type of iterator
     * IterType = Row, Col, Diag
     * Grid = constness of mm::basic_matrix
     */
    template<typename T, std::size_t Rows, std::size_t Cols, int IterType, template <typename, std::size_t, std::size_t> class Grid>
    class vector_iterator;
}

#define MM_ROW_ITER 0
#define MM_COL_ITER 1
#define MM_DIAG_ITER 2

template<typename T, std::size_t Rows, std::size_t Cols, int IterType, template <typename, std::size_t, std::size_t> class Grid>
class mm::vector_iterator
{
    std::size_t index; // variable index

    Grid<T, Rows, Cols>& M;

    const int position; // fixed index, negative too for diagonal iterator

public:
    template<typename U, std::size_t ORows, std::size_t OCols, class OIterType, template <typename, std::size_t, std::size_t> class OGrid>
    friend class vector_iterator;

    vector_iterator(Grid<T, Rows, Cols>& M, int position, std::size_t index = 0);

    mm::vector_iterator<T, Rows, Cols, IterType, Grid> operator++()
    {
        vector_iterator<T, Rows, Cols, IterType, Grid> it = *this;
        ++index;
        return it;
    }

    mm::vector_iterator<T, Rows, Cols, IterType, Grid> operator--()
    {
        vector_iterator<T, Rows, Cols, IterType, Grid> it = *this;
        --index;
        return it;
    }

    mm::vector_iterator<T, Rows, Cols, IterType, Grid>& operator++(int)
    {
        ++index;
        return *this;
    }

    mm::vector_iterator<T, Rows, Cols, IterType, Grid>& operator--(int)
    {
        --index;
        return *this;
    }

    bool operator==(const mm::vector_iterator<T, Rows, Cols, IterType, Grid>& other) const
    {
        return index == other.index;
    }

    bool operator!=(const mm::vector_iterator<T, Rows, Cols, IterType, Grid>& other) const
    {
        return index != other.index;
    }

    bool ok() const
    {
        if constexpr(IterType == MM_ROW_ITER)
            return index < Cols;
        else 
            return index < Rows;
    }

    T& operator*() const;
    T& operator[](std::size_t);
};

/* Row Iterators */

namespace mm {

    template<typename T, std::size_t Rows, std::size_t Cols>
    using row_iterator = vector_iterator<T, Rows, Cols, MM_ROW_ITER, mm::basic_matrix>;

    template<typename T, std::size_t Rows, std::size_t Cols>
    using col_iterator = vector_iterator<T, Rows, Cols, MM_COL_ITER, mm::basic_matrix>;

    template<typename T, std::size_t Rows, std::size_t Cols>
    using const_row_iterator = vector_iterator<T, Rows, Cols, MM_ROW_ITER, const mm::basic_matrix>;

    template<typename T, std::size_t Rows, std::size_t Cols>
    using const_col_iterator = vector_iterator<T, Rows, Cols, MM_COL_ITER, const mm::basic_matrix>;

    template<typename T, std::size_t N>
    using diag_iterator = vector_iterator<T, N, N, MM_DIAG_ITER, mm::basic_matrix>;

    template<typename T, std::size_t N>
    using const_diag_iterator = vector_iterator<T, N, N, MM_DIAG_ITER, const mm::basic_matrix>;
}


/*
 * Matrix class
 */

template<typename T, std::size_t Rows, std::size_t Cols>
class mm::basic_matrix {
public:
    using type = T;

    template<typename U, std::size_t ORows, std::size_t OCols>
    friend class mm::basic_matrix;

    template<typename U, std::size_t ORows, std::size_t OCols>
    friend class mm::vector_iterator;

    static constexpr std::size_t rows = Rows;
    static constexpr std::size_t cols = Cols;

    basic_matrix();

    // from initializer_list
    basic_matrix(std::initializer_list<std::initializer_list<T>> l);

    // copyable and movable
    basic_matrix(const basic_matrix<T, Rows, Cols>& other);
    basic_matrix(basic_matrix<T, Rows, Cols>&& other);

    // copy from another matrix
    template<std::size_t ORows, std::size_t OCols>
    basic_matrix(const basic_matrix<T, ORows, OCols>& other);

    // access data, basic definition
    virtual T& at(std::size_t row, std::size_t col);
    virtual const T& at(std::size_t row, std::size_t col) const;

    // allows to access a matrix M at row j col k with M[j][k]
    auto operator[](std::size_t index);

    virtual auto row_begin(std::size_t index)
    {
        return mm::row_iterator<T, Rows, Cols>(*this, static_cast<int>(index));
    }

    void swap_rows(std::size_t x, std::size_t y);
    void swap_cols(std::size_t x, std::size_t y);

    // mathematical operations
    // TODO, simply switch iteration mode
    //virtual basic_matrix<T, Cols, Rows> transposed() const;
    //inline basic_matrix<T, Cols, Rows> td() const { return transposed(); }


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
mm::basic_matrix<T, Rows, Cols>::basic_matrix(
    const mm::basic_matrix<T, Rows, Cols>& other
) : data(other.data) {}

template<typename T, std::size_t Rows, std::size_t Cols>
mm::basic_matrix<T, Rows, Cols>::basic_matrix(
    mm::basic_matrix<T, Rows, Cols>&& other
) : data(std::forward<decltype(other.data)>(other.data)) {}

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


/* member functions */

template<typename T, std::size_t Rows, std::size_t Cols>
T& mm::basic_matrix<T, Rows, Cols>::at(std::size_t row, std::size_t col) {
    assert(row < Rows); // "out of row bound"
    assert(col < Cols); // "out of column bound"

    return data[row * Cols + col];
}

template<typename T, std::size_t Rows, std::size_t Cols>
const T& mm::basic_matrix<T, Rows, Cols>::at(std::size_t row, std::size_t col) const {
    assert(row < Rows); // "out of row bound"
    assert(col < Cols); // "out of column bound"

    return data[row * Cols + col];
}

template<typename T, std::size_t Rows, std::size_t Cols>
auto mm::basic_matrix<T, Rows, Cols>::operator[](std::size_t index) {
    if constexpr (is_row_vec() || is_col_vec()) {
        return data.at(index);
    } else {
        return row_begin(index);

        /*return row_vec<T, Rows>(
            data.cbegin() + (index * Cols),
            data.cbegin() + ((index + 1) * Cols) + 1
        );*/
    }
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

    for (unsigned row = 0; row < rows; row++)
        std::swap(this->at(row, x), this->at(row, y));
}

/*template<typename T, std::size_t M, std::size_t N>
mm::basic_matrix<T, N, M> mm::basic_matrix<T, M, N>::transposed() const {
    mm::basic_matrix<T, N, M> result;

    for (unsigned row = 0; row < M; row++)
        for (unsigned col = 0; col < N; col++)
            result.at(col, row) = this->at(row, col);

    return result;
}*/


/* operator overloading */
template<typename T, std::size_t Rows, std::size_t Cols>
mm::basic_matrix<T, Rows, Cols> operator+(
    const mm::basic_matrix<T, Rows, Cols>& a,
    const mm::basic_matrix<T, Rows, Cols>& b
) {
    mm::basic_matrix<T, Rows, Cols> result;

    for (unsigned row = 0; row < Rows; row++)
        for (unsigned col = 0; col < Cols; col++)
            result.at(row, col) = a.at(row, col) + b.at(row, col);
    
    return result;
}

template<typename T, std::size_t Rows, std::size_t Cols>
mm::basic_matrix<T, Rows, Cols> operator*(
    const mm::basic_matrix<T, Rows, Cols>& m,
    const T& scalar
) {
    mm::basic_matrix<T, Rows, Cols> result;
    for (unsigned row = 0; row < Rows; row++)
        for (unsigned col = 0; col < Cols; col++)
            result.at(row, col) = m.at(row, col) * scalar;

    return result;
}

template<typename T, std::size_t Rows, std::size_t Cols>
mm::basic_matrix<T, Rows, Cols> operator*(
    const T& scalar,
    const mm::basic_matrix<T, Rows, Cols>& m
) {
    return m * scalar;
}

template<typename T, std::size_t M, std::size_t P1, std::size_t P2, std::size_t N>
mm::basic_matrix<T, M, N> operator*(
    const mm::basic_matrix<T, M, P1>& a,
    const mm::basic_matrix<T, P2, N>& b
) {
    static_assert(P1 == P2, "invalid matrix multiplication");
    mm::basic_matrix<T, M, N>  result;

    // TODO: use a more efficient algorithm
    for (unsigned row = 0; row < M; row++)
        for (unsigned col = 0; col < N; col++)
            for (unsigned k = 0; k < P1; k++)
                result.at(row, col) = a.at(row, k) * b.at(k, col);

    return result;
}


template<typename T, std::size_t Rows, std::size_t Cols>
mm::basic_matrix<T, Rows, Cols> operator-(
    const mm::basic_matrix<T, Rows, Cols>& a,
    const mm::basic_matrix<T, Rows, Cols>& b
) {
    return a + (static_cast<T>(-1) * b);
}

template<typename T, std::size_t Rows, std::size_t Cols, unsigned NumW = 3>
std::ostream& operator<<(std::ostream& os, const mm::basic_matrix<T, Rows, Cols>& m) {
    for (unsigned row = 0; row < Rows; row++) {
        os << "[ ";
        for (unsigned col = 0; col < (Cols -1); col++) {
            os << std::setw(NumW) << m.at(row, col) << ", ";
        }
        os << std::setw(NumW) << m.at(row, (Cols -1)) << " ]\n";
    }

    return os;
}



/* 
 * derivated classes 
 */

/* row vector specialization */
template<typename T, std::size_t Rows>
class mm::row_vec : public mm::basic_matrix<T, Rows, 1> {
public:
    using mm::basic_matrix<T, Rows, 1>::basic_matrix;
};

/* column vector specialization */
template<typename T, std::size_t Cols>
class mm::col_vec : public mm::basic_matrix<T, 1, Cols> {
public:
    using mm::basic_matrix<T, 1, Cols>::basic_matrix;
};

/* general specialization (alias) */
template<typename T, std::size_t Rows, std::size_t Cols>
class mm::matrix : public mm::basic_matrix<T, Rows, Cols> {
public:
    using mm::basic_matrix<T, Rows, Cols>::basic_matrix;
};

/* 
 * transposed matrix format 
 */

template<typename T, std::size_t Rows, std::size_t Cols>
class mm::transposed_matrix : public mm::basic_matrix<T, Rows, Cols>
{
public:
    using mm::basic_matrix<T, Rows, Cols>::basic_matrix;

    virtual T& at(std::size_t row, std::size_t col) override
    {
        return mm::basic_matrix<T, Rows, Cols>::at(col, row);
    }

    virtual const T& at(std::size_t row, std::size_t col) const override
    {
        return mm::basic_matrix<T, Rows, Cols>::at(col, row);
    }

    // allows to access a matrix M at row j col k with M[j][k]
    virtual auto row_begin(std::size_t index) override
    {
        return mm::col_iterator<T, Rows, Cols>(*this, static_cast<int>(index));
    }
};

/* square matrix specialization */
template<typename T, std::size_t N>
class mm::square_matrix : public mm::basic_matrix<T, N, N> {
public:
    using mm::basic_matrix<T, N, N>::basic_matrix;

    /// in place transpose
    //void transpose();  
    //inline void t() { transpose(); }

    T trace();
    inline T tr() { return trace(); }

    /// in place inverse
    // TODO, det != 0
    // TODO, use gauss jordan for invertible ones
    void invert();


    // TODO, downcast to K-diagonal, user defined cast
    /*template<int K>
    operator mm::diagonal_matrix<T, N, K>() const
    {
        // it's always possible to do it bidirectionally, 
        // without loosing information
        return reinterpret_cast<mm::diagonal_matrix<T, N, K>>(*this);
    }*/

    // get the identity of size N
    static inline constexpr square_matrix<T, N> identity() {
        square_matrix<T, N> i;
        for (unsigned row = 0; row < N; row++)
            for (unsigned col = 0; col < N; col++)
                i.at(row, col) = (row == col) ? 1 : 0;

        return i;
    }
};

/* 
 * K-diagonal square matrix format 
 * K is bounded between ]-N, N[
 */

template<typename T, std::size_t N, int K>
class mm::diagonal_matrix : public mm::square_matrix<T, N>
{
public:
    using mm::square_matrix<T, N>::square_matrix;

    // TODO, redefine at, operator[]
    // TODO, matrix multiplication
};

/*template<typename T, std::size_t N>
void mm::square_matrix<T, N>::transpose() {
    for (unsigned row = 0; row < N; row++)
        for (unsigned col = 0; col < row; col++)
            std::swap(this->at(row, col), this->at(col, row));
}*/

template<typename T, std::size_t N>
T mm::square_matrix<T, N>::trace() {

    T sum = 0;
    for (mm::diag_iterator<T, N> it(*this, 0); it.ok(); ++it)
        sum += *it;

    return sum;
}

/* Iterators implementation */

template<typename T, std::size_t Rows, std::size_t Cols, int IterType, template <typename, std::size_t, std::size_t> class Grid>
mm::vector_iterator<T, Rows, Cols, IterType, Grid>::vector_iterator(Grid<T, Rows, Cols>& _M, int pos, std::size_t i)
    : index(i), M(_M), position(pos)
{
    if constexpr (IterType == MM_ROW_ITER) {
        assert(pos < Cols);
    } else if constexpr (IterType == MM_COL_ITER) {
        assert(pos < Rows);
    } else if constexpr (IterType == MM_DIAG_ITER) {
        assert(abs(pos) < Rows);
    }
}

template<typename T, std::size_t Rows, std::size_t Cols, int IterType, template <typename, std::size_t, std::size_t> class Grid>
T& mm::vector_iterator<T, Rows, Cols, IterType, Grid>::operator*() const
{
    if constexpr (IterType == MM_ROW_ITER)
        return M.data[position * Cols + index];
    else if constexpr (IterType == MM_COL_ITER)
        return M.data[index * Cols + position];
    else if constexpr (IterType == MM_DIAG_ITER)
        return (position > 0) ?
            M.data[(index + position) * Cols + index] :
            M.data[index * Cols + (index - position)];
}

template<typename T, std::size_t Rows, std::size_t Cols, int IterType, template <typename, std::size_t, std::size_t> class Grid>
T& mm::vector_iterator<T, Rows, Cols, IterType, Grid>::operator[](std::size_t i)
{
    if constexpr (IterType == MM_ROW_ITER)
        return M.data[position * Cols + i];
    else if constexpr (IterType == MM_COL_ITER)
        return M.data[i * Cols + position];
    else if constexpr (IterType == MM_DIAG_ITER)
        return (position > 0) ?
            M.data[(i + position) * Cols + i] :
            M.data[i * Cols + (i - position)];
}

