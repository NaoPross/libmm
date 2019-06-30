#pragma once


/* 

template<typename T, std::size_t Rows, std::size_t Cols>
mm::basic_matrix<T, Rows, Cols>::basic_matrix(
    const mm::basic_matrix<T, Rows, Cols>& other
) : data(other.data) {}

template<typename T, std::size_t Rows, std::size_t Cols>
mm::basic_matrix<T, Rows, Cols>::basic_matrix(
    mm::basic_matrix<T, Rows, Cols>&& other
) : data(std::forward<decltype(other.data)>(other.data)) {}*/


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

/*template<typename T, std::size_t M, std::size_t N>
mm::basic_matrix<T, N, M> mm::basic_matrix<T, M, N>::transposed() const {
    mm::basic_matrix<T, N, M> result;

    for (unsigned row = 0; row < M; row++)
        for (unsigned col = 0; col < N; col++)
            result.at(col, row) = this->at(row, col);

    return result;
}*/
