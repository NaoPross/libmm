/* mmmatrix.hpp
 * Part of Mathematical library built (ab)using Modern C++ 17 abstractions.
 *
 * This library is not intended to be _performant_, it does not contain
 * hand written SMID / SSE / AVX optimizations. It is instead an example
 * of highly abstracted code, where matrices can contain any data type.
 *
 * As a challenge, the matrix data structure has been built on a container
 * of static capacity. But if a dynamic base container is needed, the code
 * should be easily modifiable to add further abstraction, by templating
 * the container and possibly the allocator.
 *
 * Naoki Pross <naopross@thearcway.org>
 * 2018 ~ 2019
 */
#pragma once

#include <cmath>
#include <array>

namespace mm {
    template<typename T, std::size_t Rows, std::size_t Cols>
    class basic_matrix;
}

template<typename T, std::size_t Rows, std::size_t Cols>
class mm::basic_matrix {
public:
    using type = T;

    static constexpr std::size_t rows = Rows;
    static constexpr std::size_t cols = Cols;

    basic_matrix() {}

    template<std::size_t Rows_, std::size_t Cols_>
    basic_matrix(const basic_matrix<T, Rows_, Cols_>& other);

    const T& at(std::size_t row, std::size_t col);

private:
    std::array<T, (Rows * Cols)> data;
};



template<typename T, std::size_t Rows, std::size_t Cols>
template<std::size_t ORows, std::size_t OCols>
mm::basic_matrix<T, Rows, Cols>::basic_matrix(const basic_matrix<T, ORows, OCols>& other) {
    static_assert((ORows <= Rows),
        "cannot copy a taller matrix into a smaller one"
    );

    static_assert((OCols <= Cols),
        "cannot copy a larger matrix into a smaller one"
    );

    std::copy(std::begin(other.data), std::end(other.data), data.begin());
}
