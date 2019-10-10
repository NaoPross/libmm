#pragma once

namespace mm {
    namespace algorithm {
        // does nothing
        struct visit
        {
            template<typename Matrix>
            void operator()(Matrix& m) {}
        };

        struct transpose
        {
            /// does not work with non-square matrices
            template<typename Matrix>
            void operator()(Matrix& m) {
                static_assert(Matrix::rows == Matrix::cols);
                // naiive impl
                for (index r = 0; r < m.rows / 2; r++)
                    for (index c = 0; c < m.cols; c++)
                        if (c != r)
                            std::swap(m.at(r, c), m.at(c, r));
            }
        };

        /// algorithm aliases
        using tr = transpose;
    }

    /// namespace alias
    namespace alg = algorithm;
}