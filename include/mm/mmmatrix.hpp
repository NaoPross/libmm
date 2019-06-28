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

    // TODO, not sure it's a good idea
    //template<typename T, std::size_t Rows, std::size_t Cols>
    //class transposed_matrix;

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

    template<typename T, std::size_t N>
    class diagonal_matrix;

    /* 
     * Iterators 
     */

    template<typename T, std::size_t Rows, std::size_t Cols>
    class vector_iterator;

    template<typename T, std::size_t N>
    class diag_iterator;

    template<typename T, std::size_t Rows, std::size_t Cols>
    class const_vector_iterator;

    template<typename T, std::size_t N>
    class const_diag_iterator;
}

/* Non-const Iterators */

template<typename T, std::size_t Rows, std::size_t Cols>
class mm::vector_iterator
{
    std::size_t index; // variable index

    mm::basic_matrix<T, Rows, Cols>& M;

    const std::size_t position; // fixed index
    const bool direction; // true = row, false = column

public:
    template<typename U, std::size_t ORows, std::size_t OCols>
    friend class vector_iterator;

    vector_iterator(mm::basic_matrix<T, Rows, Cols>& M, std::size_t position, bool direction);

    mm::vector_iterator<T, Rows, Cols> operator++()
    {
        vector_iterator<T, Rows, Cols> it = *this;
        ++index;
        return it;
    }

    mm::vector_iterator<T, Rows, Cols> operator--()
    {
        vector_iterator<T, Rows, Cols> it = *this;
        --index;
        return it;
    }

    mm::vector_iterator<T, Rows, Cols>& operator++(int)
    {
        ++index;
        return *this;
    }

    mm::vector_iterator<T, Rows, Cols>& operator--(int)
    {
        --index;
        return *this;
    }

    bool operator==(const mm::vector_iterator<T, Rows, Cols>& other) const
    {
        return index == other.index;
    }

    bool operator=!(const mm::vector_iterator<T, Rows, Cols>& other) const
    {
        return index != other.index;
    }

    T& operator*() const;
    T& operator[](std::size_t);
};

template<typename T, std::size_t N>
class diag_iterator
{
    std::size_t index; // variable index

    mm::square_matrix<T, N>& M;

    const int position; // fixed diagonal index

public:
    template<typename U, std::size_t ON>
    friend class diag_iterator;

    diag_iterator(mm::square_matrix<T, N>& M, std::size_t position, bool direction);

    mm::diag_iterator<T, N> operator++()
    {
        diag_iterator<T, N> it = *this;
        ++index;
        return it;
    }

    mm::diag_iterator<T, N> operator--()
    {
        diag_iterator<T, N> it = *this;
        --index;
        return it;
    }

    mm::diag_iterator<T, N>& operator++(int)
    {
        ++index;
        return *this;
    }

    mm::diag_iterator<T, N>& operator--(int)
    {
        --index;
        return *this;
    }

    bool operator==(const mm::diag_iterator<T, N>& other) const
    {
        return index == other.index;
    }

    bool operator=!(const mm::diag_iterator<T, N>& other) const
    {
        return index != other.index;
    }

    T& operator*() const;
};

/* Const Iterators */

template<typename T, std::size_t Rows, std::size_t Cols>
class mm::const_vector_iterator
{
    std::size_t index; // variable index

    const mm::basic_matrix<T, Rows, Cols>& M;

    const std::size_t position; // fixed index
    const bool direction; // true = row, false = column

public:
    const_vector_iterator(mm::basic_matrix<T, Rows, Cols>& M, std::size_t position, bool direction);

    mm::const_vector_iterator<T, Rows, Cols> operator++()
    {
        vector_iterator<T, Rows, Cols> it = *this;
        ++index;
        return it;
    }

    mm::const_vector_iterator<T, Rows, Cols> operator--()
    {
        vector_iterator<T, Rows, Cols> it = *this;
        --index;
        return it;
    }

    mm::const_vector_iterator<T, Rows, Cols>& operator++(int)
    {
        ++index;
        return *this;
    }

    mm::const_vector_iterator<T, Rows, Cols>& operator--(int)
    {
        --index;
        return *this;
    }

    bool operator==(const mm::const_vector_iterator<T, Rows, Cols>& other) const
    {
        return index == other;
    }

    bool operator=!(const mm::const_vector_iterator<T, Rows, Cols>& other) const
    {
        return index != other;
    }

    const T& operator*() const;
    const T& operator[](std::size_t) const;
};

template<typename T>
class const_diag_iterator
{
    std::size_t index; // variable index

    const mm::square_matrix<T, N>& M;

    const int position; // fixed diagonal index

public:
    template<typename U, std::size_t ON>
    friend class const_diag_iterator;

    const_diag_iterator(const mm::square_matrix<T, N>& M, std::size_t position, bool direction);

    mm::const_diag_iterator<T, N> operator++()
    {
        const_diag_iterator<T, N> it = *this;
        ++index;
        return it;
    }

    mm::const_diag_iterator<T, N> operator--()
    {
        const_diag_iterator<T, N> it = *this;
        --index;
        return it;
    }

    mm::const_diag_iterator<T, N>& operator++(int)
    {
        ++index;
        return *this;
    }

    mm::const_diag_iterator<T, N>& operator--(int)
    {
        --index;
        return *this;
    }

    bool operator==(const mm::const_diag_iterator<T, N>& other) const
    {
        return index == other.index;
    }

    bool operator=!(const mm::const_diag_iterator<T, N>& other) const
    {
        return index != other.index;
    }

    const T& operator*() const;
};

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

    // access data
    virtual T& at(std::size_t row, std::size_t col);
    virtual const T& at(std::size_t row, std::size_t col) const;

    // allows to access a matrix M at row j col k with M[j][k]
    virtual auto operator[](std::size_t index);

    void swap_rows(std::size_t x, std::size_t y);
    void swap_cols(std::size_t x, std::size_t y);

    // mathematical operations
    // TODO, simply switch iteration mode
    virtual basic_matrix<T, Cols, Rows> transposed() const;
    inline basic_matrix<T, Cols, Rows> td() const { return transposed(); }


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
        return row_vec<T, Rows>(
            data.cbegin() + (index * Cols),
            data.cbegin() + ((index + 1) * Cols) + 1
        );
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

template<typename T, std::size_t M, std::size_t N>
mm::basic_matrix<T, N, M> mm::basic_matrix<T, M, N>::transposed() const {
    mm::basic_matrix<T, N, M> result;

    for (unsigned row = 0; row < M; row++)
        for (unsigned col = 0; col < N; col++)
            result.at(col, row) = this->at(row, col);

    return result;
}


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
 * TODO: write this class, or put a bool flag into the original one
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
    virtual auto operator[](std::size_t index) override
    {
        // TODO, return other direction iterator
    }
}

/* square matrix specialization */
template<typename T, std::size_t N>
class mm::square_matrix : public mm::basic_matrix<T, N, N> {
public:
    using mm::basic_matrix<T, N, N>::basic_matrix;

    /// in place transpose
    void transpose();  
    inline void t() { transpose(); }

    T trace();
    inline T tr() { return trace(); }

    /// in place inverse
    // TODO, det != 0
    // TODO, use gauss jordan for invertible ones
    void invert();


    // TODO, downcast to K-diagonal, user defined cast
    template<int K>
    operator mm::diagonal_matrix<T, N, K>() const
    {
        // it's always possible to do it bidirectionally, 
        // without loosing information
        return dynamic_cast<mm::diagonal_matrix<T, N, K>>(*this);
    }

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
class mm::diagonal_matrix : public mm::square_matrix
{
public:
    using mm::square_matrix<T, N>::square_matrix;

    // TODO, redefine at, operator[]
    // TODO, matrix multiplication
};

template<typename T, std::size_t N>
void mm::square_matrix<T, N>::transpose() {
    for (unsigned row = 0; row < N; row++)
        for (unsigned col = 0; col < row; col++)
            std::swap(this->at(row, col), this->at(col, row));
}

template<typename T, std::size_t N>
T mm::square_matrix<T, N>::trace() {
    T sum = 0;
    for (unsigned i = 0; i < N; i++)
        sum += this->at(i, i);

    return sum;
}

/* Iterators implementations */

template<typename T, std::size_t Rows, std::size_t Cols>
mm::vector_iterator<T, Rows, Cols>::vector_iterator(mm::basic_matrix<T, Rows, Cols>& _M, std::size_t pos, bool dir)
    index(0), M(_M), position(pos), direction(dir)
{
    assert((dir && pos < Cols) || (!dir && pos < Rows))
}

template<typename T, std::size_t Rows, std::size_t Cols>
T& mm::vector_iterator<T, Rows, Cols>::operator*() const
{
    return (direction) ?
            M.data[position * Cols + index] :
            M.data[index * Cols + position];
}

template<typename T, std::size_t Rows, std::size_t Cols>
T& mm::vector_iterator<T, Rows, Cols>::operator[](std::size_t i)
{
    return (direction) ?
            M.data[position * Cols + i] :
            M.data[i * Cols + position];
}

template<typename T, std::size_t N>
mm::diag_iterator<T, N>::diag_iterator(mm::square_matrix<T, N>& _M, int pos)
    index(0), M(_M), position(pos)
{
    assert(abs(pos) < N) // pos bounded between ]-N, N[ 
}

template<typename T, std::size_t N>
T& mm::diag_iterator<T, N>::operator*() const
{
    return (k > 0) ?
        M.data[(index + position) * Cols + index] :
        M.data[index * Cols + (index - position)];
}


template<typename T, std::size_t Rows, std::size_t Cols>
mm::const_vector_iterator<T, Rows, Cols>::const_vector_iterator(const mm::basic_matrix<T, Rows, Cols>& _M, std::size_t pos, bool dir)
    index(0), M(_M), position(pos), direction(dir)
{
    assert((dir && pos < Cols) || (!dir && pos < Rows))
}

template<typename T, std::size_t Rows, std::size_t Cols>
const T& mm::const_vector_iterator<T, Rows, Cols>::operator*() const
{
    return (direction) ?
            M.data[position * Cols + index] :
            M.data[index * Cols + position];
}

template<typename T, std::size_t Rows, std::size_t Cols>
const T& mm::const_vector_iterator<T, Rows, Cols>::operator[](std::size_t i) const
{
    return (direction) ?
            M.data[position * Cols + i] :
            M.data[i * Cols + position];
}

template<typename T, std::size_t N>
mm::const_diag_iterator<T, N>::const_diag_iterator(const mm::square_matrix<T, N>& _M, int pos)
    index(0), M(_M), position(pos)
{
    assert(abs(pos) < N) // pos bounded between ]-N, N[ 
}

template<typename T, std::size_t N>
T& mm::const_diag_iterator<T, N>::operator*() const
{
    return (k > 0) ?
        M.data[(index + position) * Cols + index] :
        M.data[index * Cols + (index - position)];
}



