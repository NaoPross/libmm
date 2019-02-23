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
#include <cstring>
#include <cassert>
#include <initializer_list>

namespace mm {
    template<typename T, std::size_t Rows, std::size_t Cols>
    class basic_matrix;

    template<typename T, std::size_t Rows, std::size_t Cols>
    class matrix;

    template<typename T, std::size_t N>
    class square_matrix;

    // template<typename T, std::size_t N>
    // class diag_matrix;

    template<typename T, std::size_t Rows>
    class row_vec;

    template<typename T, std::size_t Cols>
    class col_vec;
}

template<typename T, std::size_t Rows, std::size_t Cols>
class mm::basic_matrix {
public:
    using type = T;

    static constexpr std::size_t rows = Rows;
    static constexpr std::size_t cols = Cols;

    // from initializer_list
    basic_matrix(std::initializer_list<std::initializer_list<T>> l);

    // copyable and movable
    basic_matrix(const basic_matrix<T, Rows, Cols>& other);
    basic_matrix(basic_matrix<T, Rows, Cols>&& other);

    // copy from another matrix
    template<std::size_t ORows, std::size_t OCols>
    basic_matrix(const basic_matrix<T, ORows, OCols>& other);

    // copy or move from 2D array
    basic_matrix(const T (& values)[Rows][Cols]);
    basic_matrix(T (&& values)[Rows][Cols]);

    // access data
    T& at(std::size_t row, std::size_t col);
    const T& at(std::size_t row, std::size_t col) const;
    auto&& operator[](std::size_t index);

    void swap_rows(std::size_t x, std::size_t y);
    void swap_cols(std::size_t x, std::size_t y);

    // mathematical operations
    basic_matrix<T, Cols, Rows> transposed() const;
    inline basic_matrix<T, Cols, Rows> trd() const { return transposed(); }

    // bool is_invertible() const;
    // basic_matrix<T, Rows, Cols> inverse() const;


    /// downcast to square matrix
    inline constexpr bool is_square() const { return (Rows == Cols); }
    inline constexpr square_matrix<T, Rows> to_square() const {
        static_assert(is_square());
        return static_cast<square_matrix<T, Rows>>(*this);
    }


    /// downcast to row_vector
    inline constexpr bool is_row_vec() const { return (Cols == 1); }
    inline constexpr row_vec<T, Rows> to_row_vec() const {
        static_assert(is_row_vec());
        return static_cast<row_vec<T, Rows>>(*this);
    }

    /// downcast to col_vector
    inline constexpr bool is_col_vec() const { return (Rows == 1); }
    inline constexpr col_vec<T, Cols> to_col_vec() const {
        static_assert(is_col_vec());
        return static_cast<col_vec<T, Cols>>(*this);
    }

private:
    T data[Rows * Cols] = {};
};


template<typename T, std::size_t Rows, std::size_t Cols>
mm::basic_matrix<T, Rows, Cols>::basic_matrix(
    std::initializer_list<std::initializer_list<T>> l
) {
    assert(l.size() == Rows);

    auto row_it = l.begin();
    for (unsigned row = 0; (row < Rows) && (row_it != l.end()); row++) {
        assert((*row_it).size() == Cols);

        auto col_it = (*row_it).begin();
        for (unsigned col = 0; (col < Cols) && (col_it != (*row_it).end()); col++) {
            this->at(row, col) = *col_it;
            ++col_it;
        }
        ++row_it;
    }
}

template<typename T, std::size_t Rows, std::size_t Cols>
mm::basic_matrix<T, Rows, Cols>::basic_matrix(
    const mm::basic_matrix<T, Rows, Cols>& other
) {
    std::memcpy(&data, &other.data, sizeof(data));
}

template<typename T, std::size_t Rows, std::size_t Cols>
mm::basic_matrix<T, Rows, Cols>::basic_matrix(
    mm::basic_matrix<T, Rows, Cols>&& other
) {
    data = other.data;
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

    std::memcpy(&data, &other.data, sizeof(data));
}

template<typename T, std::size_t Rows, std::size_t Cols>
mm::basic_matrix<T, Rows, Cols>::basic_matrix(const T (& values)[Rows][Cols]) {
    std::memcpy(&data, &values, sizeof(data));
}

template<typename T, std::size_t Rows, std::size_t Cols>
mm::basic_matrix<T, Rows, Cols>::basic_matrix(T (&& values)[Rows][Cols]) {
    data = values;
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
auto&& mm::basic_matrix<T, Rows, Cols>::operator[](std::size_t index) {
    if constexpr (is_row_vec()) {
        static_assert(index < Rows);
        return data[index];
    } else if constexpr (is_col_vec()) {
        static_assert(index < Cols);
        return data[index];
    }

    // TODO: fix
    // return row_vec<T, Rows>(std::move(data[index]));
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
            result.at(row, col) = this->at(row, col);

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
            result.at(row, col) = a.at(row, col) + a.at(row, col);
    
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

template<typename T, std::size_t M, std::size_t P, std::size_t N>
mm::basic_matrix<T, M, N> operator*(
    const mm::basic_matrix<T, M, P>& a,
    const mm::basic_matrix<T, P, N>& b
) {
    mm::basic_matrix<T, M, N>  result;

    // TODO: use a more efficient algorithm
    for (unsigned row = 0; row < M; row++)
        for (unsigned col = 0; col < N; col++)
            for (int k = 0; k < P; k++)
                result.at(row, col) = a.at(row, k) * b.at(k, col);

    return result;
}


template<typename T, std::size_t Rows, std::size_t Cols>
mm::basic_matrix<T, Rows, Cols> operator-(
    const mm::basic_matrix<T, Rows, Cols>& a,
    const mm::basic_matrix<T, Rows, Cols>& b
) {
    return a + static_cast<T>(-1) * b;
}

template<typename T, std::size_t Rows, std::size_t Cols>
std::ostream& operator<<(std::ostream& os, const mm::basic_matrix<T, Rows, Cols>& m) {
    for (unsigned row = 0; row < Rows; row++) {
        os << "[ ";
        for (unsigned col = 0; col < (Cols -1); col++) {
            os << m.at(row, col) << ", ";
        }
        os << m.at(row, (Cols -1)) << " ]\n";
    }

    return os;
}


/* square matrix specializaiton */

template<typename T, std::size_t N>
class mm::square_matrix : public mm::basic_matrix<T, N, N> {
public:
    /// in place transpose
    void transpose();  
    inline void tr() { transpose(); }

    /// in place inverse
    void invert();
};


template<typename T, std::size_t N>
void mm::square_matrix<T, N>::transpose() {
    for (unsigned row = 0; row < N; row++)
        for (unsigned col = 0; col < row; col++)
            std::swap(this->at(row, col), this->at(col, row));
}


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

template<typename T, std::size_t Rows, std::size_t Cols>
class mm::matrix : public mm::basic_matrix<T, Rows, Cols> {
public:
    using mm::basic_matrix<T, Rows, Cols>::basic_matrix;
};
