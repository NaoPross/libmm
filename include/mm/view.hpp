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
        std::tuple<decltype(alg::visit), Algs...> visitors;

        explicit mutate(Matrix& m) : matrix(m), visitors(alg::visit) {}

        template<typename ...OAlgs, typename Alg>
        explicit mutate(Matrix& m, std::tuple<OAlgs...>&& t, Alg&& v)
            : matrix(m), visitors(std::tuple_cat(t, std::make_tuple(v))) {}

        ~mutate() {
            visit();
        }

        void visit() {
            std::apply([&](auto&&... v) {
                (v(matrix),...);
            }, visitors);
        }

        operator Matrix() {
            return std::move(matrix);
        }
    };

    template<typename Matrix, typename Alg>
    clone<Matrix> operator|(clone<Matrix>&& cl, Alg&& v) {
        /// must be callable with argument Matrix&
        static_assert(std::is_invocable<Alg, Matrix&>::value);

        /// apply alg operator
        v(cl.matrix);

        /// forward to next alg
        return std::move(cl);
    }

    template<typename Matrix, typename ...Algs, typename Alg>
    mutate<Matrix, Algs..., Alg> operator|(mutate<Matrix, Algs...>&& mut, Alg&& v) {
        /// must be callable with argument Matrix&
        static_assert(std::is_invocable<Alg, Matrix&>::value);

        /// append alg to the visitors tuple
        return mutate<Matrix, Algs..., Alg>(mut.matrix, std::move(mut.visitors), v);
    }

    /// operator| between matrix and alg lambda
    /// defaults to wrapping matrix with a mutable
    template<typename Matrix, typename Alg>
    mutate<Matrix, Alg> operator|(Matrix& m, Alg& v) {
        return mutate(m) | std::move(v);
    }

    /// operator| between matrix and alg rvalue object
    /// defaults to wrapping matrix with a mutable
    template<typename Matrix, typename Alg>
    mutate<Matrix, Alg> operator|(Matrix& m, Alg&& v) {
        return m | static_cast<Alg&>(v);
    }
}
