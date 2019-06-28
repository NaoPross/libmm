#pragma once

namespace mm {

    template<typename T>
    class diag_component;

    template<typename T, std::size_t N>
    class multi_diag_matrix;
}

/*
 * Optimized case of square matrix
 * It's a matrix only composed by a diagonal 
 */

template<class T>
class mm::diag_component
{
public:
    virtual int dimension() const = 0;
};

template<class T, std::size_t N>
class mm::diag_vector
{
public:

    // TODO, define constructor

    virtual int dimension() const override
    {
        return N;
    }

private:
    std::array<T, N - ((Diag < 0) ? -Diag : Diag)> vector;
};

template<typename T, std::size_t N>
class mm::multi_diag_matrix {
public:
    using type = T;

    template<typename U, std::size_t N>
    friend class mm::multi_diag_matrix;

    multi_diag_matrix() : shared_zero(0) {}
    ~multi_diag_matrix();

    // copyable and movable
    multi_diag_matrix(const multi_diag_matrix<T, N>& other);
    multi_diag_matrix(multi_diag_matrix<T, N>&& other);

    // copy from another matrix
    template<std::size_t N>
    multi_diag_matrix(const multi_diag_matrix<T, N>& other);

    // standard access data
    T& at(std::size_t row, std::size_t col);
    const T& at(std::size_t row, std::size_t col) const;

    // allows to access a matrix M at row j col k with M[j][k]
    auto operator[](std::size_t index);

    // swap two diagonals
    void swap_diags(std::size_t k, std::size_t l);

    // diagonal construction or substitution
    template<int Diag, int K = N - ((Diag < 0) ? -Diag : Diag)>
    void put_diag(const mm::diag_vector<T, K>& diag)
    {
        //static_assert((Diag <= -N) || (Diag >= N),
        static_assert(K < 1,
            "Diagonal number must be bounded between ]-N,N[")

        auto exist = diagonals.find(Diag);

        if (exist != diagonals.end())
            // copy
            *exists = diag;
        else
            // create and copy
            diagonals.insert(new mm::diag_vector<T, K>(diag));
    }

    // mathematical operations
    virtual multi_diag_matrix<T, N> transposed() const;
    inline multi_diag_matrix<T, N> td() const { return transposed(); }

    // multiplication rhs and lhs
    // TODO, need super class matrix abstraction and auto return type

    // A * M, TODO abstraction virtual method
    template <std::size_t Rows>
    basic_matrix<Rows, N> rhs_mult(const mm::basic_matrix<T, Rows, N>& A) const;

    // M * A, TODO abstraction virtual method
    template <std::size_t Cols>
    basic_matrix<N, Cols> lhs_mult(const mm::basic_matrix<T, N, Cols>& A) const;

protected:
    template<typename ConstIterator>
    multi_diag_matrix(ConstIterator begin, ConstIterator end);

private:
    // return an arbitrary zero in non-const mode
    T shared_zero;

    // ordered set of diagonals
    std::unordered_map<int, mm::diag_component<T>*> diagonals;
};

template<typename T, std::size_t N>
T& mm::multi_diag_matrix<T, N>::at(std::size_t row, std::size_t col) {
    assert(row < N); // "out of row bound"
    assert(col < N); // "out of column bound"

    const int k = row - col; 
    auto diag = diagonals.find(k);
    const int line = (k > 0) ? col : row;

    return (diag == diagonals.end()) ? (shared_zero = 0) : (*diag)[line];
}

template<typename T, std::size_t N>
const T& mm::multi_diag_matrix<T, N>::at(std::size_t row, std::size_t col) const {
    assert(row < N); // "out of row bound"
    assert(col < N); // "out of column bound"

    const int k = row - col;
    auto diag = diagonals.find(k);
    const int line = (k > 0) ? col : row;

    return (diag == diagonals.end()) ? 0 : (*diag)[line];
}

template<typename T, std::size_t N>
auto mm::multi_diag_matrix<T, N>::operator[](std::size_t index) {
    assert(index < N)
    
    // TODO, single row mapping
}

template <typename T, std::size_t N, std::size_t Rows>
mm::basic_matrix<Rows, N> mm::multi_diag_matrix<T, N>::rhs_mult(const mm::basic_matrix<T, Rows, N>& A) const
{
    // TODO
}

template <typename T, std::size_t N, std::size_t Cols>
mm::basic_matrix<N, Cols> mm::multi_diag_matrix<T, N>::lhs_mult(const mm::basic_matrix<T, N, Cols>& A) const
{
    mm::basic_matrix<N, Cols> out;


}


