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

    /* specialization of basic_matrix for Rows == Cols */
    template<typename T, std::size_t N>
    class square_interface;

    /*
     * Most general type of matrix iterator
     * IterType = Row, Column, Diagonal
     * Grid = constness of mm::basic_matrix
     */
    template<typename T, std::size_t Rows, std::size_t Cols, int IterType, class Grid>
    class vector_iterator;

    /*
     * Access methods interface
     */

    template<typename T, std::size_t Rows, std::size_t Cols, bool Regular>
    class access;

    /* specialisations */

    template<typename T, std::size_t Rows, std::size_t Cols>
    class matrix; // simple matrix format

    template<typename T, std::size_t Rows, std::size_t Cols>
    class t_matrix; // transposed matrix format

    /* specialization of basic_matrx for Cols = 1 */
    template<typename T, std::size_t Rows>
    class row_vec;

    /* specialization of basic_matrx for Rows = 1 */
    template<typename T, std::size_t Cols>
    class col_vec; // transposed version of row_vec

    template<typename T, std::size_t N>
    class square_matrix;

    template<typename T, std::size_t N>
    class t_square_matrix;

    /* specialisation of a square_matrix for a sub-diagonal composed matrix */
    template<typename T, std::size_t N, std::size_t K = 0>
    class diag_matrix;

    template<typename T, std::size_t N, std::size_t K = 0>
    class t_diag_matrix;
}

// TODO, short term solution
#define MM_ROW_ITER 0
#define MM_COL_ITER 1
#define MM_DIAG_ITER 2

template<typename T, std::size_t Rows, std::size_t Cols, int IterType, class Grid>
class mm::vector_iterator
{
    std::size_t index; // variable index

    Grid& M;

    const int position; // fixed index, negative too for diagonal iterator

public:
    template<typename U, std::size_t ORows, std::size_t OCols, class OIterType, class OGrid>
    friend class vector_iterator;

    vector_iterator(Grid& M, int position, std::size_t index = 0);

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
    using row_iterator = vector_iterator<T, Rows, Cols, MM_ROW_ITER, mm::basic_matrix<T, Rows, Cols>>;

    template<typename T, std::size_t Rows, std::size_t Cols>
    using col_iterator = vector_iterator<T, Rows, Cols, MM_COL_ITER, mm::basic_matrix<T, Rows, Cols>>;

    template<typename T, std::size_t Rows, std::size_t Cols>
    using const_row_iterator = vector_iterator<T, Rows, Cols, MM_ROW_ITER, std::add_const<mm::basic_matrix<T, Rows, Cols>>>;

    template<typename T, std::size_t Rows, std::size_t Cols>
    using const_col_iterator = vector_iterator<T, Rows, Cols, MM_COL_ITER, std::add_const<mm::basic_matrix<T, Rows, Cols>>>;

    template<typename T, std::size_t N>
    using diag_iterator = vector_iterator<T, N, N, MM_DIAG_ITER, mm::basic_matrix<T, N, N>>;

    template<typename T, std::size_t N>
    using const_diag_iterator = vector_iterator<T, N, N, MM_DIAG_ITER, std::add_const<mm::basic_matrix<T, N, N>>>;
}

/*
 * Accessors
 */
template<typename T, std::size_t Rows, std::size_t Cols, bool Regular>
class mm::access
{
public:
    
    //access(mm::basic_matrix<T, Rows, Cols>& ref) : M(ref) {}

    T& at(std::size_t row, std::size_t col);
    const T& at(std::size_t row, std::size_t col) const;

    auto operator[](std::size_t index);
    auto operator[](std::size_t index) const;

//private:
//    mm::basic_matrix<T, Rows, Cols>& M;
protected:
    std::array<T, Rows * Cols> data;
};

/*
 * Square interface
 */

template<typename T, std::size_t N>
class mm::square_interface {
public:

    //square_interface(mm:basic_matrix<T, N, N>& _M) : M(_M) {}

    T trace();
    inline T tr() { return trace(); }

    // TODO, determinant

    /// in place inverse
    // TODO, det != 0
    // TODO, use gauss jordan for invertible ones
    //void invert();, TODO, section algorithm

//private:
//    mm:basic_matrix<T, N, N>& M; // one information more than mm::matrix !
protected:
    std::array<T, N * N> data;
};


/*
 * Matrix class, no access methods
 */

template<typename T, std::size_t Rows, std::size_t Cols>
class mm::basic_matrix {
public:
    using type = T;

    template<typename U, std::size_t ORows, std::size_t OCols>
    friend class mm::basic_matrix;

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
    //virtual T& at(std::size_t row, std::size_t col);
    //virtual const T& at(std::size_t row, std::size_t col) const; 

    // allows to access a matrix M at row j col k with M[j][k]
    //auto operator[](std::size_t index);

    void swap_rows(std::size_t x, std::size_t y);
    void swap_cols(std::size_t x, std::size_t y);

    // mathematical operations
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

//private:
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

/*template<typename T, std::size_t Rows, std::size_t Cols>
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
        return mm::row_iterator<T, Rows, Cols>(*this, static_cast<int>(index));

        //return row_vec<T, Rows>(
         //   data.cbegin() + (index * Cols),
          //  data.cbegin() + ((index + 1) * Cols) + 1
        );
    }
}*/


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


/* TODO, operator overloading */
/*template<typename T, std::size_t Rows, std::size_t Cols>
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
}*/

/* 
 * derivated classes 
 */

// simple format matrix
template<typename T, std::size_t Rows, std::size_t Cols>
class mm::matrix : public mm::basic_matrix<T, Rows, Cols>, virtual public mm::access<T, Rows, Cols, true>
{
public:
    using mm::basic_matrix<T, Rows, Cols>::basic_matrix;
};

// transposed matrix
template<typename T, std::size_t Rows, std::size_t Cols>
class mm::t_matrix : public mm::basic_matrix<T, Rows, Cols>, virtual public mm::access<T, Rows, Cols, false>
{
public:
    using mm::basic_matrix<T, Rows, Cols>::basic_matrix;
};

/* row vector specialization */
template<typename T, std::size_t Rows>
class mm::row_vec : public mm::matrix<T, Rows, 1> {
public:
    using mm::matrix<T, Rows, 1>::matrix;

    // TODO, begin, end
};

/* column vector specialization */
template<typename T, std::size_t Cols>
class mm::col_vec : public mm::t_matrix<T, 1, Cols> {
public:
    using mm::t_matrix<T, 1, Cols>::t_matrix;

    // TODO, begin, end
};

/* square matrix specialization */
template<typename T, std::size_t N>
class mm::square_matrix : public mm::matrix<T, N, N> , virtual public mm::square_interface<T, N> {
public:
    using mm::matrix<T, N, N>::matrix;

    // get the identity of size N
    static inline constexpr square_matrix<T, N> identity() {
        mm::square_matrix<T, N> i;
        for (unsigned row = 0; row < N; row++)
            for (unsigned col = 0; col < N; col++)
                i.at(row, col) = (row == col) ? 1 : 0;

        return i;
    }
};

template<typename T, std::size_t N>
class mm::t_square_matrix : virtual public mm::t_matrix<T, N, N>, virtual public mm::square_interface<T, N> {
public:
    using mm::t_matrix<T, N, N>::t_matrix;
};

/* 
 * K-diagonal square matrix format 
 * K is bounded between ]-N, N[
 */

/*template<typename T, std::size_t N, std::size_t K>
class mm::diagonal_matrix : public mm::square_matrix<T, N>
{
public:
    using mm::square_matrix<T, N>::square_matrix;

    // TODO, redefine at, operator[]
    // TODO, matrix multiplication
};*/

/*template<typename T, std::size_t N>
void mm::square_matrix<T, N>::transpose() {
    for (unsigned row = 0; row < N; row++)
        for (unsigned col = 0; col < row; col++)
            std::swap(this->at(row, col), this->at(col, row));
}*/

/* Iterators implementation */

template<typename T, std::size_t Rows, std::size_t Cols, int IterType, class Grid>
mm::vector_iterator<T, Rows, Cols, IterType, Grid>::vector_iterator(Grid& _M, int pos, std::size_t i)
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

template<typename T, std::size_t Rows, std::size_t Cols, int IterType, class Grid>
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

template<typename T, std::size_t Rows, std::size_t Cols, int IterType, class Grid>
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

/*
 * Accessors implementation
 */

template<typename T, std::size_t Rows, std::size_t Cols, bool Regular>
T& mm::access<T, Rows, Cols, Regular>::at(std::size_t row, std::size_t col)
{
    if constexpr (Regular)
        return data[row * Cols + col];
    else
        return data[col * Cols + row]; // transpose
}

template<typename T, std::size_t Rows, std::size_t Cols, bool Regular>
const T& mm::access<T, Rows, Cols, Regular>::at(std::size_t row, std::size_t col) const
{
    if constexpr (Regular)
        return data[row * Cols + col];
    else
        return data[col * Cols + row]; // transpose
}

template<typename T, std::size_t Rows, std::size_t Cols, bool Regular>
auto mm::access<T, Rows, Cols, Regular>::operator[](std::size_t index)
{
    if constexpr (this->is_row_vec() || this->is_col_vec())
        return data.at(index);
    else if (Regular)
        return mm::row_iterator<T, Rows, Cols>(*this, static_cast<int>(index));
    else
        return mm::col_iterator<T, Rows, Cols>(*this, static_cast<int>(index));
}

template<typename T, std::size_t Rows, std::size_t Cols, bool Regular>
auto mm::access<T, Rows, Cols, Regular>::operator[](std::size_t index) const
{
    if constexpr (this->is_row_vec() || this->is_col_vec())
        return data.at(index);
    else if (Regular)
        return mm::const_row_iterator<T, Rows, Cols>(*this, static_cast<int>(index));
    else
        return mm::const_col_iterator<T, Rows, Cols>(*this, static_cast<int>(index));
}

/* Square interface implementation */

template<typename T, std::size_t N>
T mm::square_interface<T, N>::trace() 
{
    T sum = 0;
    for (mm::diag_iterator<T, N> it(*this, 0); it.ok(); ++it)
        sum += *it;

    return sum;
}

