#pragma once

#include "mm/mmmatrix.hpp"

namespace mm {
    namespace view  {
        namespace details {
            struct viewer_base {};
        }

        template<typename Matrix, enable_if_matrix_t<Matrix>>
        struct row
        {
            using type = typename Matrix::type;
            using cols = typename Matrix::cols;

            Matrix& matrix;
            index i;

            explicit row(index r, Matrix& m) : i(r), matrix(m) {}

            type& at(index j) {
                return matrix.at(i, j);
            }
        };

        template<typename Matrix, enable_if_matrix_t<Matrix>>
        struct col
        {
            using type = typename Matrix::type;
            using rows = typename Matrix::rows;

            Matrix& matrix;
            index j;

            explicit col(index c, Matrix& m) : j(c), matrix(m) {}

            type& at(index i) {
                return matrix.at(i, j);
            }
        };

        template<typename Matrix,
            std::size_t SubRows,
            std::size_t SubCols,
            enable_if_matrix_t<Matrix>>
        struct submatrix 
        {
            using type = typename Matrix::type;
            static constexpr std::size_t rows = SubRows;
            static constexpr std::size_t cols = SubCols;

            Matrix& matrix;
            index i, j;

            explicit submatrix(index row, index col, Matrix& m) : i(row), j(col), matrix(m) {}
        };
    }

    namespace algorithm {
        /// does nothing
        constexpr auto visit = [](auto& m) {};

        /* transpose */

        // TODO: fix impl
        constexpr auto transpose = [](auto& m) {
            using Matrix = typename std::decay<decltype(m)>::type;
            static_assert(Matrix::rows == Matrix::cols);

            // naiive impl
            for (index r = 0; r < Matrix::rows / 2; r++)
                for (index c = 0; c < Matrix::cols; c++)
                    if (c != r)
                        std::swap(m.at(r, c), m.at(c, r));
        };

        constexpr auto& tr = transpose;

        /* inverse matrix */
        template<typename Matrix, enable_if_matrix_t<Matrix>>
        struct inverse_of
        {
            Matrix matrix;

            inverse_of(const Matrix& m) : matrix(m) {}

            void operator()(Matrix& m) {
                // TODO
            }

            operator Matrix() {
                std::move(matrix);
            }
        };

        constexpr auto invert = [](auto& m) { 
            return inverse_of(m);
        };

        /* row reduction */
        template<typename Matrix, typename _en = enable_if_matrix_t<Matrix>>
        struct row_reduction_of
        {
            using T = typename Matrix::type;

            Matrix& matrix;

            row_reduction_of(Matrix& m) : matrix(m) {}

            void compute() {
                index pivot_col = 0;
                for (index r = 0; r < Matrix::rows; r++) {
                    auto rv = view::row(r, matrix);

                    // divide row
                    for (index c = 0; c < Matrix::cols; c++) {
                        // rv.at(c) /= pivot;
                    }

                    // next row - row[0] * prev row
                }
            }

            operator Matrix() {
                compute();
                return std::move(matrix);
            }

            void operator()(Matrix& m) {
                compute();
                // m = std::move(matrix);
            }
        };

        template<typename Matrix>
        using rref_of = row_reduction_of<Matrix>;

        constexpr auto rref = [](auto& m) {
            return row_reduction_of(m);
        };
    }


    /// namespace alias
    namespace alg = algorithm;
}
