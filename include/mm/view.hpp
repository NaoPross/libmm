#pragma once

#include "mm/mmmatrix.hpp"
#include "mm/algorithm.hpp"

#include <tuple>
#include <type_traits>


namespace mm {
    template<typename Matrix>
    struct clone
    {
        Matrix matrix;

        explicit clone(Matrix&& m) : matrix(m) {}
        explicit clone(const Matrix& m) : matrix(m) {}

        operator Matrix() {
            return std::move(matrix);
        }
    };

    template<typename Matrix, typename ...Algs>
    struct mutate
    {
        Matrix& matrix;
        std::tuple<Algs...> visitors;

        explicit mutate(Matrix& m) : matrix(m) {}

        template<typename ...OAlgs, typename Alg>
        explicit mutate(Matrix& m, std::tuple<OAlgs...>&& t, Alg&& v)
            : matrix(m)
        {
            /// append the new operator
            visitors = std::tuple_cat(t, std::make_tuple(v));
        }

        ~mutate() {
            visit();
        }

        void visit() {
            std::apply([this](auto&&... v) {
                (v(matrix),...);
            }, visitors);
        }

        operator Matrix() {
            return std::move(matrix);
        }
    };

    template<typename Matrix, typename Alg>
    clone<Matrix> operator|(clone<Matrix>&& cl, Alg&& v) {
        /// apply alg operator
        v(cl.matrix);
        /// forward to next alg
        return clone<Matrix>(std::move(cl));
    }

    template<typename Matrix, typename ...Algs, typename Alg>
    mutate<Matrix, Algs..., Alg> operator|(mutate<Matrix, Algs...>&& mut, Alg&& v) {
        /// append alg to the visitors tuple
        return mutate<Matrix, Algs..., Alg>(
            mut.matrix,
            std::move(mut.visitors),
            v
        );
    }

    template<typename Matrix, typename Alg>
    mutate<Matrix, Alg> operator|(Matrix& m, Alg&& v) {
        return mutate(m) | std::move(v);
    }
}
