#pragma once

#include "debug.hpp"

namespace mm::iter {

    template<typename T, std::size_t Rows, std::size_t Cols, class IterType, class Grid>
    class vector_iterator;

    template<typename T, std::size_t Rows, std::size_t Cols, class Grid>
    class basic_iterator;

    template<typename T, std::size_t N, class Grid>
    class diag_iterator;
}

template<typename T, std::size_t Rows, std::size_t Cols, class IterType, class Grid>
class mm::iter::vector_iterator
{
public:

    template<typename U, std::size_t R, std::size_t C, class G>
    friend class mm::iter::basic_iterator;

    template<typename U, std::size_t N, class G>
    friend class mm::iter::diag_iterator;

    vector_iterator(Grid& _M, std::size_t pos, std::size_t i = 0)
        : M(_M), position(pos), index(i) {}

//#ifdef MM_IMPLICIT_CONVERSION_ITERATOR
    operator T&()
    {
        //npdebug("Calling +")
        return *(*this);
    }
//#endif

    IterType operator++()
    {
        IterType it = cpy();
        ++index;
        return it;
    }

    IterType operator--()
    {
        IterType it = cpy();
        --index;
        return it;
    }

    IterType& operator++(int)
    {
        ++index;
        return ref();
    }

    IterType& operator--(int)
    {
        --index;
        return ref();
    }

    bool operator==(const IterType& other) const
    {
        return index == other.index;
    }

    bool operator!=(const IterType& other) const
    {
        return index != other.index;
    }

    bool ok() const
    {
        return index < size();
    }

    virtual std::size_t size() const = 0;
    
    virtual T& operator*() = 0;
    virtual T& operator[](std::size_t) = 0;

    virtual T& operator[](std::size_t) const = 0;

    IterType begin()
    {
        return IterType(M, position, 0);
    }

    virtual IterType end() = 0;
    
protected:

    Grid& M; // grid mapping

    const std::size_t position; // fixed index, negative too for diagonal iterator
    std::size_t index; // variable index

    virtual IterType& ref() = 0;
    virtual IterType cpy() = 0;
};


/*
 * Scalar product 
 */

template<typename T, 
    std::size_t R1, std::size_t C1, 
    std::size_t R2, std::size_t C2, 
    class IterType1, class IterType2,
    class Grid1, class Grid2>
typename std::remove_const<T>::type operator*(const mm::iter::vector_iterator<T, R1, C1, IterType1, Grid1>& v,
            const mm::iter::vector_iterator<T, R2, C2, IterType2, Grid2>& w)
{
    typename std::remove_const<T>::type out(0);
    const std::size_t N = std::min(v.size(), w.size());

    for(unsigned i = 0; i < N; ++i)
        out += v[i] * w[i];

    return out;
}

template<typename T, std::size_t Rows, std::size_t Cols, class Grid>
class mm::iter::basic_iterator : public mm::iter::vector_iterator<T, Rows, Cols, mm::iter::basic_iterator<T, Rows, Cols, Grid>, Grid>
{
    bool direction;

    virtual mm::iter::basic_iterator<T, Rows, Cols, Grid>& ref() override
    {
        return *this;
    }

    virtual mm::iter::basic_iterator<T, Rows, Cols, Grid> cpy() override
    {
        return *this;
    }

public:

    basic_iterator(Grid& A, std::size_t pos, std::size_t _index = 0, bool dir = true)
        : mm::iter::vector_iterator<T, Rows, Cols, mm::iter::basic_iterator<T, Rows, Cols, Grid>, Grid>
            (A, pos, _index), direction(dir)
    {
        //npdebug("Position: ", pos, ", Rows: ", Rows, " Cols: ", Cols, ", Direction: ", dir)

        if (direction)
           assert(pos < Rows);
        else
           assert(pos < Cols);
    }

    virtual std::size_t size() const
    {
        return (direction) ? Cols : Rows;
    }

    
    virtual T& operator*() override
    {
        return (direction) ?
            this->M.data[this->position * Cols + this->index] :
            this->M.data[this->index * Cols + this->position];

    }

    virtual T& operator[](std::size_t i) override
    {
        return (direction) ?
            this->M.data[this->position * Cols + i] :
            this->M.data[i * Cols + this->position];
    }

    virtual T& operator[](std::size_t i) const override
    {
        return (direction) ?
            this->M.data[this->position * Cols + i] :
            this->M.data[i * Cols + this->position];
    }

    virtual mm::iter::basic_iterator<T, Rows, Cols, Grid> end()
    {
        return mm::iter::basic_iterator<T, Rows, Cols, Grid>(this->M, this->position, 
                    (direction) ? Cols : Rows);
    }
};

template<typename T, std::size_t N, class Grid>
class mm::iter::diag_iterator : public mm::iter::vector_iterator<T, N, N, mm::iter::diag_iterator<T, N, Grid>, Grid>
{
    bool sign;

    virtual mm::iter::diag_iterator<T, N, Grid>& ref() override
    {
        return *this;
    }

    virtual mm::iter::diag_iterator<T, N, Grid> cpy() override
    {
        return *this;
    }

public:

    diag_iterator(Grid& A, signed long int pos, std::size_t _index = 0)
        : mm::iter::vector_iterator<T, N, N, mm::iter::diag_iterator<T, N, Grid>, Grid>
            (A, static_cast<std::size_t>(labs(pos)), _index), sign(pos >= 0)
    {
        assert(this->position < N);
    }

    virtual std::size_t size() const
    {
        return N - this->position;
    }

    virtual T& operator*() override
    {
        return (sign) ?
            this->M.data[(this->index - this->position) * N + this->index] :
            this->M.data[this->index * N + (this->index + this->position)];
    }

    virtual T& operator[](std::size_t i) override
    {
        return (sign) ?
            this->M.data[(i - this->position) * N + i] :
            this->M.data[i * N + (i + this->position)];
    }

    virtual T& operator[](std::size_t i) const override
    {
        return (sign) ?
            this->M.data[(i - this->position) * N + i] :
            this->M.data[i * N + (i + this->position)];
    }


    virtual mm::iter::diag_iterator<T, N, Grid> end()
    {
        return mm::iter::diag_iterator<T, N, Grid>(this->M, this->position, N);
    }
};


