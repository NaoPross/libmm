#pragma once

template<typename T, std::size_t Rows, std::size_t Cols, bool Regular>
class mm::access : virtual public mm::basic_matrix<T, Rows, Cols>
{
public:
    template<typename U, std::size_t ORows, std::size_t OCols>
    friend class mm::access;

    using mm::basic_matrix<T, Rows, Cols>::basic_matrix;

    access(mm::basic_matrix<T, Rows, Cols>&& m) 
        : basic_matrix<T, Rows, Cols>(m) {}

    virtual T& at(std::size_t row, std::size_t col) = 0;
    virtual const T& at(std::size_t row, std::size_t col) const = 0;

    static constexpr std::size_t rows = Regular ? Rows : Cols;
    static constexpr std::size_t cols = Regular ? Cols : Rows;
};

template<typename T, std::size_t Rows, std::size_t Cols>
class mm::row_access : public mm::access<T, Rows, Cols, true>
{
public:
    template<typename U, std::size_t ORows, std::size_t OCols>
    friend class mm::row_access;

    using mm::access<T, Rows, Cols, true>::access;

    row_access(mm::basic_matrix<T, Rows, Cols>&& m) 
        : access<T, Rows, Cols, true>(m) {}

    virtual T& at(std::size_t row, std::size_t col) override;
    virtual const T& at(std::size_t row, std::size_t col) const override;

    row_iterator<T, Rows, Cols> operator[](std::size_t index);
    const_row_iterator<T, Rows, Cols> operator[](std::size_t index) const;
};

template<typename T, std::size_t Rows, std::size_t Cols>
class mm::col_access : public mm::access<T, Rows, Cols, false>
{
public:
    template<typename U, std::size_t ORows, std::size_t OCols>
    friend class mm::access;

    using mm::access<T, Rows, Cols, false>::access;

    col_access(mm::basic_matrix<T, Rows, Cols>&& m) 
        : access<T, Rows, Cols, false>(m) {}

    virtual T& at(std::size_t row, std::size_t col) override;
    virtual const T& at(std::size_t row, std::size_t col) const override;

    col_iterator<T, Rows, Cols> operator[](std::size_t index);
    const_col_iterator<T, Rows, Cols> operator[](std::size_t index) const;
};

/*
 * Square interface
 */

template<typename T, std::size_t N>
class mm::square_interface : virtual public mm::basic_matrix<T, N, N>
{
public:

    template<typename U, std::size_t ON>
    friend class mm::square_interface;

    using mm::basic_matrix<T, N, N>::basic_matrix;

    T trace();
    inline T tr() { return trace(); }

    mm::diag_iterator<T, N> diagonal_iter(int index = 0)
    {
        return mm::diag_iterator<T, N>(*this, index);
    }

    mm::diag_iterator<T, N> diagonal_iter(int index = 0) const
    {
        return mm::const_diag_iterator<T, N>(*this, index);
    }

    // TODO, determinant

    /// in place inverse
    // TODO, det != 0
    // TODO, use gauss jordan for invertible ones
    //void invert();, TODO, section algorithm
};

/* 
 * derivated classes 
 */

// simple format matrix
template<typename T, std::size_t Rows, std::size_t Cols>
class mm::matrix : public mm::row_access<T, Rows, Cols>
{
public:
    using mm::row_access<T, Rows, Cols>::row_access;

    matrix(mm::t_matrix<T, Rows, Cols>&& m)
        : row_access<T, Rows, Cols>(std::move<mm::basic_matrix<T, Rows, Cols>>(m))
    {
                
    }

    operator mm::t_matrix<T, Rows, Cols>()
    {
        return mm::t_matrix<T, Rows, Cols>(*this); 
    }

    mm::t_matrix<T, Rows, Cols>& t()
    {
        return static_cast<mm::t_matrix<T, Rows, Cols>>(*this);
    }

    const mm::t_matrix<T, Rows, Cols>& t() const
    {
        return static_cast<mm::t_matrix<T, Rows, Cols>>(*this);
    }
};

template<typename T, std::size_t Rows, std::size_t Cols>
class mm::t_matrix : public mm::col_access<T, Rows, Cols>
{
public:
    using mm::col_access<T, Rows, Cols>::col_access;

    t_matrix(mm::matrix<T, Rows, Cols>&& m)
        : col_access<T, Rows, Cols>(std::move<mm::basic_matrix<T, Rows, Cols>>(m))
    {
                
    }

    operator mm::matrix<T, Rows, Cols>()
    {
        return mm::matrix<T, Rows, Cols>(*this); 
    }

    mm::matrix<T, Rows, Cols>& t()
    {
        return static_cast<mm::matrix<T, Rows, Cols>>(*this);
    }

    const mm::matrix<T, Rows, Cols>& t() const
    {
        return static_cast<mm::matrix<T, Rows, Cols>>(*this);
    }
};

// transposed matrix
/*template<typename T, std::size_t Rows, std::size_t Cols>
class mm::t_matrix : public mm::basic_matrix<T, Rows, Cols>, virtual public mm::access<T, Rows, Cols, false>
{
public:
    using mm::basic_matrix<T, Rows, Cols>::basic_matrix;
};*/

/* row vector specialization */
template<typename T, std::size_t Rows>
class mm::row_vec : public mm::matrix<T, Rows, 1> {
public:
    using mm::matrix<T, Rows, 1>::matrix;

    // TODO, begin, end
};

/* column vector specialization */
template<typename T, std::size_t Cols>
class mm::col_vec : public mm::matrix<T, 1, Cols> {
public:
    using mm::matrix<T, 1, Cols>::matrix;

    // TODO, begin, end
};

/* square matrix specialization */
template<typename T, std::size_t N>
class mm::square_matrix : public mm::matrix<T, N, N>, public mm::square_interface<T, N> {
public:
    using mm::matrix<T, N, N>::matrix;

    // get the identity of size N
    static inline constexpr square_matrix<T, N> identity() {
        mm::square_matrix<T, N> i;
        for (unsigned row = 0; row < N; row++)
            for (unsigned col = 0; col < N; col++)
                i.at(row, col) = (row == col) ? 1 : 0;

        return i;
    }
};

template<typename T, std::size_t N>
class mm::t_square_matrix : public mm::t_matrix<T, N, N>, public mm::square_interface<T, N> {
public:
    using mm::t_matrix<T, N, N>::t_matrix;
};

/* 
 * K-diagonal square matrix format 
 * K is bounded between ]-N, N[
 */

/*template<typename T, std::size_t N, std::size_t K>
class mm::diagonal_matrix : public mm::square_matrix<T, N>
{
public:
    using mm::square_matrix<T, N>::square_matrix;

    // TODO, redefine at, operator[]
    // TODO, matrix multiplication
};*/

/*
 * Accessors implementation
 */

template<typename T, std::size_t Rows, std::size_t Cols>
T& mm::row_access<T, Rows, Cols>::at(std::size_t row, std::size_t col)
{
    return mm::basic_matrix<T, Rows, Cols>::data[row * Cols + col];
}

template<typename T, std::size_t Rows, std::size_t Cols>
const T& mm::row_access<T, Rows, Cols>::at(std::size_t row, std::size_t col) const
{
    return mm::basic_matrix<T, Rows, Cols>::data[row * Cols + col];
}

template<typename T, std::size_t Rows, std::size_t Cols>
T& mm::col_access<T, Rows, Cols>::at(std::size_t row, std::size_t col)
{
    return mm::basic_matrix<T, Rows, Cols>::data[col * Cols + row]; // transpose
}

template<typename T, std::size_t Rows, std::size_t Cols>
const T& mm::col_access<T, Rows, Cols>::at(std::size_t row, std::size_t col) const
{
    return mm::basic_matrix<T, Rows, Cols>::data[col * Cols + row]; // transpose
}

template<typename T, std::size_t Rows, std::size_t Cols>
mm::row_iterator<T, Rows, Cols> mm::row_access<T, Rows, Cols>::operator[](std::size_t index)
{
    return mm::row_iterator<T, Rows, Cols>(*this, static_cast<int>(index));
}

template<typename T, std::size_t Rows, std::size_t Cols>
mm::const_row_iterator<T, Rows, Cols> mm::row_access<T, Rows, Cols>::operator[](std::size_t index) const
{
    return mm::const_row_iterator<T, Rows, Cols>(*this, static_cast<int>(index));
}

template<typename T, std::size_t Rows, std::size_t Cols>
mm::col_iterator<T, Rows, Cols> mm::col_access<T, Rows, Cols>::operator[](std::size_t index)
{
    return mm::col_iterator<T, Rows, Cols>(*this, static_cast<int>(index));
}

template<typename T, std::size_t Rows, std::size_t Cols>
mm::const_col_iterator<T, Rows, Cols> mm::col_access<T, Rows, Cols>::operator[](std::size_t index) const
{
    return mm::const_col_iterator<T, Rows, Cols>(*this, static_cast<int>(index));
}

/* Square interface implementation */

template<typename T, std::size_t N>
T mm::square_interface<T, N>::trace() 
{
    T sum = 0;
    for (const auto& x : diagonal_iter())
        sum += x;

    return sum;
}


