#pragma once

namespace mm {
    namespace view  {
        template<typename Matrix>
        struct row
        {
            using T = typename Matrix::type;

            Matrix& matrix;
            index i;

            explicit row(index r, Matrix& m) : i(r), matrix(m) {}

            T& at(index j) {
                return matrix.at(i, j);
            }
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
        template<typename Matrix>
        struct inverse_of
        {
            Matrix matrix;

            inverse_of(const Matrix& m) : matrix(m) {}

            void operator()(Matrix& m) {
                // TODO
            }
        };

        constexpr auto invert = [](auto& m) { 
            return inverse_of(m);
        };

        /* row reduction */
        template<typename Matrix>
        struct row_reduction_of
        {
            using T = typename Matrix::type;

            Matrix matrix;

            row_reduction_of(const Matrix& m) : matrix(m) {}

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