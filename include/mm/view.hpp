#pragma once

#include <mmmatrix.hpp>


namespace mm::alg {

	template <
		template<typename, std::size_t, std::size_t> typename Matrix,
		typename T, std::size_t Rows, std::size_t Cols
	>
	struct visitor
	{
		using type = T;

		// copy constructible
		visitor(const visitor& other) = default;

		T& operator()(const Matrix<T, Rows, Cols>& m, index row, index col) {
			return m.at(row, col);
		}

		const T& operator()(const Matrix<T, Rows, Cols>& m, index row, index col) {
			return operator()(m, row, col);
		}
	};

	template <
		template<typename, std::size_t, std::size_t> typename Matrix,
		typename T, std::size_t Rows, std::size_t Cols
	>
	struct transpose : public visitor<Matrix, T, Rows, Cols>
	{
		T& operator()(const Matrix<T, Rows, Cols> m, index row, index col) {
			// assert(col < Rows)
			// assert(row < Cols)
			return m.at(col, row);
		}
	};
}

namespace mm {
	template <
		template<typename, std::size_t, std::size_t> typename Matrix,
		typename T, std::size_t Rows, std::size_t Cols
	>
	struct view
	{
		Matrix<T, Rows, Cols>& m;
		// std::stack<std::unique_ptr<alg::visitor>> visitors;
		std::unique_ptr<alg::visitor> visitor;

		T& at(index row, index col) {
			return visitor(m, row, col);
		}

		view& operator|=(const alg::visitor& other) {
			// visitors.push(std::move(std::make_unique<alg::visitor>(other)));
			visitor = std::make_unique<alg::visistor>(other);
		}
	};

	view operator|(const view& left, const alg::visitor& right) {
		return left |= right;
	}
}
