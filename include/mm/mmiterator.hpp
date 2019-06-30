#pragma once

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

    operator T&()
    {
        return *(*this);
    }

    IterType operator++()
    {
        IterType it = *this;
        ++index;
        return it;
    }

    IterType operator--()
    {
        IterType it = *this;
        --index;
        return it;
    }

    IterType& operator++(int)
    {
        ++index;
        return *this;
    }

    IterType& operator--(int)
    {
        --index;
        return *this;
    }

    bool operator==(const IterType& other) const
    {
        return index == other.index;
    }

    bool operator!=(const IterType& other) const
    {
        return index != other.index;
    }

    virtual bool ok() const = 0;
    
    virtual T& operator*() = 0;
    virtual T& operator[](std::size_t) = 0;

    IterType begin()
    {
        return IterType(M, position, 0);
    }

    virtual IterType end() = 0;
    
    /*
     * Scalar product 
     */

    template<std::size_t P>
    T operator*(const mm::iter::vector_iterator<T, Rows, P, IterType, Grid>& v)
    {
        T out(0);

        for (unsigned k(0); k < Rows; ++k)
            out += (*this)[k] * v[k];

        return out;
    }

    template<std::size_t P>
    T operator*(const mm::iter::vector_iterator<T, P, Cols, IterType, Grid>& v)
    {
        T out(0);

        for (unsigned k(0); k < Cols; ++k)
            out += (*this)[k] * v[k];

        return out;
    }

protected:

    Grid& M; // grid mapping

    const std::size_t position; // fixed index, negative too for diagonal iterator
    std::size_t index; // variable index
};

template<typename T, std::size_t Rows, std::size_t Cols, class Grid>
class mm::iter::basic_iterator : public mm::iter::vector_iterator<T, Rows, Cols, mm::iter::basic_iterator<T, Rows, Cols, Grid>, Grid>
{
    bool direction;

public:

    basic_iterator(Grid& A, std::size_t pos, std::size_t _index = 0, bool dir = true)
        : mm::iter::vector_iterator<T, Rows, Cols, mm::iter::basic_iterator<T, Rows, Cols, Grid>, Grid>
            (A, pos, _index), direction(dir)
    {
        if (direction)
           assert(pos < Rows);
        else
           assert(pos < Cols);
    }

    virtual bool ok() const override
    {
        return (direction) ? this->index < Cols : this->index < Rows;
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

public:

    diag_iterator(Grid& A, signed long pos, std::size_t _index = 0)
        : mm::iter::vector_iterator<T, N, N, mm::iter::diag_iterator<T, N, Grid>, Grid>
            (A, static_cast<std::size_t>(abs(pos)), _index), sign(pos >= 0)
    {
        assert(this->position < N);
    }

    virtual bool ok() const 
    {
        return this->index < N;
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

    virtual mm::iter::diag_iterator<T, N, Grid> end()
    {
        return mm::iter::diag_iterator<T, N, Grid>(this->M, this->position, N);
    }
};


