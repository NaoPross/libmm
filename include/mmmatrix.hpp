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

namespace mm {
    template<typename T, std::size_t Rows, std::size_t Cols>
    class basic_matrix;

    template<typename T, std::size_t Rows>
    using row_vec = basic_matrix<T, Rows, 1>;

    template<typename T, std::size_t Cols>
    using col_vec = basic_matrix<T, 1, Cols>;
}

template<typename T, std::size_t Rows, std::size_t Cols>
class mm::basic_matrix {
public:
    using type = T;

    static constexpr std::size_t rows = Rows;
    static constexpr std::size_t cols = Cols;

    basic_matrix(const basic_matrix<T, Rows, Cols>& other);
    basic_matrix(basic_matrix<T, Rows, Cols>&& other);

    template<std::size_t ORows, std::size_t OCols>
    basic_matrix(const basic_matrix<T, ORows, OCols>& other);

    // access data
    T& at(std::size_t row, std::size_t col);

    void swap_rows(std::size_t x, std::size_t y);
    void swap_cols(std::size_t x, std::size_t y);

    // mathematical operations
    basic_matrix<T, Cols, Rows> transposed();
    inline basic_matrix<T, Cols, Rows> t() { return transposed(); }

    // bool is_invertible();
    // bool invert();
    // basic_matrix<T, Rows, Cols> inverse();

    inline constexpr bool is_square() {
        return (Rows == Cols);
    }

private:
    T data[Rows][Cols] = {};
};


template<typename T, std::size_t Rows, std::size_t Cols>
mm::basic_matrix<T, Rows, Cols>::basic_matrix(const mm::basic_matrix<T, Rows, Cols>& other) {
    for (int row = 0; row < Rows; row++)
        for (int col = 0; col < Cols; col++)
            data[row][col] = other.data[row][col];
}

template<typename T, std::size_t Rows, std::size_t Cols>
mm::basic_matrix<T, Rows, Cols>::basic_matrix(mm::basic_matrix<T, Rows, Cols>&& other) {
    data = other.data;
}

template<typename T, std::size_t Rows, std::size_t Cols>
template<std::size_t ORows, std::size_t OCols>
mm::basic_matrix<T, Rows, Cols>::basic_matrix(const mm::basic_matrix<T, ORows, OCols>& other) {
    static_assert((ORows <= Rows),
        "cannot copy a taller matrix into a smaller one"
    );

    static_assert((OCols <= Cols),
        "cannot copy a larger matrix into a smaller one"
    );

    for (int row = 0; row < Rows; row++)
        for (int col = 0; col < Cols; col++)
            data[row][col] = other.data[row][col];
}


/* member functions */

template<typename T, std::size_t Rows, std::size_t Cols>
T& mm::basic_matrix<T, Rows, Cols>::at(std::size_t row, std::size_t col) {
    static_assert(row < Rows, "out of row bound");
    static_assert(col < Cols, "out of column bound");

    return data[row][col];
}

template<typename T, std::size_t Rows, std::size_t Cols>
void mm::basic_matrix<T, Rows, Cols>::swap_rows(std::size_t x, std::size_t y) {
    if (x == y)
        return;

    for (int col = 0; col < Cols; col++)
        std::swap(data[x][col], data[y][col]);
}

template<typename T, std::size_t Rows, std::size_t Cols>
void mm::basic_matrix<T, Rows, Cols>::swap_cols(std::size_t x, std::size_t y) {
    if (x == y)
        return;

    for (int row = 0; row < rows; row++)
        std::swap(data[row][x], data[row][y]);
}

template<typename T, std::size_t M, std::size_t N>
mm::basic_matrix<T, N, M> mm::basic_matrix<T, M, N>::transposed() {
    mm::basic_matrix<T, N, M> result;

    for (int row = 0; row < M; row++)
        for (int col = 0; col < N; col++)
            result.at(row, col) = at(col, row);

    return result;
}


/* operator overloading */
template<typename T, std::size_t Rows, std::size_t Cols>
mm::basic_matrix<T, Rows, Cols> operator+(
    const mm::basic_matrix<T, Rows, Cols>& a,
    const mm::basic_matrix<T, Rows, Cols>& b
) {
    mm::basic_matrix<T, Rows, Cols> result;

    for (int row = 0; row < Rows; row++)
        for (int col = 0; col < Cols; col++)
            result.at(row, col) = a.at(row, col) + a.at(row, col);
    
    return result;
}

template<typename T, std::size_t Rows, std::size_t Cols>
mm::basic_matrix<T, Rows, Cols> operator*(
    const mm::basic_matrix<T, Rows, Cols>& m,
    const T& scalar
) {
    mm::basic_matrix<T, Rows, Cols> result;
    for (int row = 0; row < Rows; row++)
        for (int col = 0; col < Cols; col++)
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
    for (int row = 0; row < M; row++)
        for (int col = 0; col < N; col++)
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
    for (int row = 0; row < Rows; row++) {
        os << "[ ";
        for (int col = 0; col < (Cols -1); col++) {
            os << m.at(row, col);
        }
        os << m.at(Rows -1, Cols -1) << " ]\n";
    }

    return os;
}
