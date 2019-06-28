#pragma once

template<typename T, std::size_t Rows, std::size_t Cols>
class mm::vector_iterator
{
    std::size_t index; // variable index

    mm::basic_matrix<T, Rows, Cols>& M;

    const std::size_t position; // fixed index
    const bool direction; // true = row, false = column

public:
    template<typename U, std::size_t ORows, std::size_t OCols>
    friend class vector_iterator;

    vector_iterator(mm::basic_matrix<T, Rows, Cols>& M, std::size_t position, bool direction);

    mm::vector_iterator<T, Rows, Cols> operator++()
    {
        vector_iterator<T, Rows, Cols> it = *this;
        ++index;
        return it;
    }

    mm::vector_iterator<T, Rows, Cols> operator--()
    {
        vector_iterator<T, Rows, Cols> it = *this;
        --index;
        return it;
    }

    mm::vector_iterator<T, Rows, Cols>& operator++(int)
    {
        ++index;
        return *this;
    }

    mm::vector_iterator<T, Rows, Cols>& operator--(int)
    {
        --index;
        return *this;
    }

    bool operator==(const mm::vector_iterator<T, Rows, Cols>& other) const
    {
        return index == other.index;
    }

    bool operator=!(const mm::vector_iterator<T, Rows, Cols>& other) const
    {
        return index != other.index;
    }

    T& operator*() const;
    T& operator[](std::size_t);
};

template<typename T, std::size_t N>
class diag_iterator
{
    std::size_t index; // variable index

    mm::square_matrix<T, N>& M;

    const int position; // fixed diagonal index

public:
    template<typename U, std::size_t ON>
    friend class diag_iterator;

    diag_iterator(mm::square_matrix<T, N>& M, std::size_t position, bool direction);

    mm::diag_iterator<T, N> operator++()
    {
        diag_iterator<T, N> it = *this;
        ++index;
        return it;
    }

    mm::diag_iterator<T, N> operator--()
    {
        diag_iterator<T, N> it = *this;
        --index;
        return it;
    }

    mm::diag_iterator<T, N>& operator++(int)
    {
        ++index;
        return *this;
    }

    mm::diag_iterator<T, N>& operator--(int)
    {
        --index;
        return *this;
    }

    bool operator==(const mm::diag_iterator<T, N>& other) const
    {
        return index == other.index;
    }

    bool operator=!(const mm::diag_iterator<T, N>& other) const
    {
        return index != other.index;
    }

    T& operator*() const;
};

/* Const Iterators */

template<typename T, std::size_t Rows, std::size_t Cols>
class mm::const_vector_iterator
{
    std::size_t index; // variable index

    const mm::basic_matrix<T, Rows, Cols>& M;

    const std::size_t position; // fixed index
    const bool direction; // true = row, false = column

public:
    const_vector_iterator(mm::basic_matrix<T, Rows, Cols>& M, std::size_t position, bool direction);

    mm::const_vector_iterator<T, Rows, Cols> operator++()
    {
        vector_iterator<T, Rows, Cols> it = *this;
        ++index;
        return it;
    }

    mm::const_vector_iterator<T, Rows, Cols> operator--()
    {
        vector_iterator<T, Rows, Cols> it = *this;
        --index;
        return it;
    }

    mm::const_vector_iterator<T, Rows, Cols>& operator++(int)
    {
        ++index;
        return *this;
    }

    mm::const_vector_iterator<T, Rows, Cols>& operator--(int)
    {
        --index;
        return *this;
    }

    bool operator==(const mm::const_vector_iterator<T, Rows, Cols>& other) const
    {
        return index == other;
    }

    bool operator=!(const mm::const_vector_iterator<T, Rows, Cols>& other) const
    {
        return index != other;
    }

    const T& operator*() const;
    const T& operator[](std::size_t) const;
};

template<typename T>
class const_diag_iterator
{
    std::size_t index; // variable index

    const mm::square_matrix<T, N>& M;

    const int position; // fixed diagonal index

public:
    template<typename U, std::size_t ON>
    friend class const_diag_iterator;

    const_diag_iterator(const mm::square_matrix<T, N>& M, std::size_t position, bool direction);

    mm::const_diag_iterator<T, N> operator++()
    {
        const_diag_iterator<T, N> it = *this;
        ++index;
        return it;
    }

    mm::const_diag_iterator<T, N> operator--()
    {
        const_diag_iterator<T, N> it = *this;
        --index;
        return it;
    }

    mm::const_diag_iterator<T, N>& operator++(int)
    {
        ++index;
        return *this;
    }

    mm::const_diag_iterator<T, N>& operator--(int)
    {
        --index;
        return *this;
    }

    bool operator==(const mm::const_diag_iterator<T, N>& other) const
    {
        return index == other.index;
    }

    bool operator=!(const mm::const_diag_iterator<T, N>& other) const
    {
        return index != other.index;
    }

    const T& operator*() const;
};

/* Iterators implementations */

template<typename T, std::size_t Rows, std::size_t Cols>
mm::vector_iterator<T, Rows, Cols>::vector_iterator(mm::basic_matrix<T, Rows, Cols>& _M, std::size_t pos, bool dir)
    index(0), M(_M), position(pos), direction(dir)
{
    assert((dir && pos < Cols) || (!dir && pos < Rows))
}

template<typename T, std::size_t Rows, std::size_t Cols>
T& mm::vector_iterator<T, Rows, Cols>::operator*() const
{
    return (direction) ?
            M.data[position * Cols + index] :
            M.data[index * Cols + position];
}

template<typename T, std::size_t Rows, std::size_t Cols>
T& mm::vector_iterator<T, Rows, Cols>::operator[](std::size_t i)
{
    return (direction) ?
            M.data[position * Cols + i] :
            M.data[i * Cols + position];
}

template<typename T, std::size_t N>
mm::diag_iterator<T, N>::diag_iterator(mm::square_matrix<T, N>& _M, int pos)
    index(0), M(_M), position(pos)
{
    assert(abs(pos) < N) // pos bounded between ]-N, N[ 
}

template<typename T, std::size_t N>
T& mm::diag_iterator<T, N>::operator*() const
{
    return (k > 0) ?
        M.data[(index + position) * Cols + index] :
        M.data[index * Cols + (index - position)];
}


template<typename T, std::size_t Rows, std::size_t Cols>
mm::const_vector_iterator<T, Rows, Cols>::const_vector_iterator(const mm::basic_matrix<T, Rows, Cols>& _M, std::size_t pos, bool dir)
    index(0), M(_M), position(pos), direction(dir)
{
    assert((dir && pos < Cols) || (!dir && pos < Rows))
}

template<typename T, std::size_t Rows, std::size_t Cols>
const T& mm::const_vector_iterator<T, Rows, Cols>::operator*() const
{
    return (direction) ?
            M.data[position * Cols + index] :
            M.data[index * Cols + position];
}

template<typename T, std::size_t Rows, std::size_t Cols>
const T& mm::const_vector_iterator<T, Rows, Cols>::operator[](std::size_t i) const
{
    return (direction) ?
            M.data[position * Cols + i] :
            M.data[i * Cols + position];
}

template<typename T, std::size_t N>
mm::const_diag_iterator<T, N>::const_diag_iterator(const mm::square_matrix<T, N>& _M, int pos)
    index(0), M(_M), position(pos)
{
    assert(abs(pos) < N) // pos bounded between ]-N, N[ 
}

template<typename T, std::size_t N>
T& mm::const_diag_iterator<T, N>::operator*() const
{
    return (k > 0) ?
        M.data[(index + position) * Cols + index] :
        M.data[index * Cols + (index - position)];
}



/* Iterators implementations */

template<typename T, std::size_t Rows, std::size_t Cols>
mm::vector_iterator<T, Rows, Cols>::vector_iterator(mm::basic_matrix<T, Rows, Cols>& _M, std::size_t pos, bool dir)
    index(0), M(_M), position(pos), direction(dir)
{
    assert((dir && pos < Cols) || (!dir && pos < Rows))
}

template<typename T, std::size_t Rows, std::size_t Cols>
T& mm::vector_iterator<T, Rows, Cols>::operator*() const
{
    return (direction) ?
            M.data[position * Cols + index] :
            M.data[index * Cols + position];
}

template<typename T, std::size_t Rows, std::size_t Cols>
T& mm::vector_iterator<T, Rows, Cols>::operator[](std::size_t i)
{
    return (direction) ?
            M.data[position * Cols + i] :
            M.data[i * Cols + position];
}

template<typename T, std::size_t N>
mm::diag_iterator<T, N>::diag_iterator(mm::square_matrix<T, N>& _M, int pos)
    index(0), M(_M), position(pos)
{
    assert(abs(pos) < N) // pos bounded between ]-N, N[ 
}

template<typename T, std::size_t N>
T& mm::diag_iterator<T, N>::operator*() const
{
    return (k > 0) ?
        M.data[(index + position) * Cols + index] :
        M.data[index * Cols + (index - position)];
}


template<typename T, std::size_t Rows, std::size_t Cols>
mm::const_vector_iterator<T, Rows, Cols>::const_vector_iterator(const mm::basic_matrix<T, Rows, Cols>& _M, std::size_t pos, bool dir)
    index(0), M(_M), position(pos), direction(dir)
{
    assert((dir && pos < Cols) || (!dir && pos < Rows))
}

template<typename T, std::size_t Rows, std::size_t Cols>
const T& mm::const_vector_iterator<T, Rows, Cols>::operator*() const
{
    return (direction) ?
            M.data[position * Cols + index] :
            M.data[index * Cols + position];
}

template<typename T, std::size_t Rows, std::size_t Cols>
const T& mm::const_vector_iterator<T, Rows, Cols>::operator[](std::size_t i) const
{
    return (direction) ?
            M.data[position * Cols + i] :
            M.data[i * Cols + position];
}

template<typename T, std::size_t N>
mm::const_diag_iterator<T, N>::const_diag_iterator(const mm::square_matrix<T, N>& _M, int pos)
    index(0), M(_M), position(pos)
{
    assert(abs(pos) < N) // pos bounded between ]-N, N[ 
}

template<typename T, std::size_t N>
T& mm::const_diag_iterator<T, N>::operator*() const
{
    return (k > 0) ?
        M.data[(index + position) * Cols + index] :
        M.data[index * Cols + (index - position)];
}


/* Iterators implementations */

template<typename T, std::size_t Rows, std::size_t Cols>
mm::vector_iterator<T, Rows, Cols>::vector_iterator(mm::basic_matrix<T, Rows, Cols>& _M, std::size_t pos, bool dir)
    index(0), M(_M), position(pos), direction(dir)
{
    assert((dir && pos < Cols) || (!dir && pos < Rows))
}

template<typename T, std::size_t Rows, std::size_t Cols>
T& mm::vector_iterator<T, Rows, Cols>::operator*() const
{
    return (direction) ?
            M.data[position * Cols + index] :
            M.data[index * Cols + position];
}

template<typename T, std::size_t Rows, std::size_t Cols>
T& mm::vector_iterator<T, Rows, Cols>::operator[](std::size_t i)
{
    return (direction) ?
            M.data[position * Cols + i] :
            M.data[i * Cols + position];
}

template<typename T, std::size_t N>
mm::diag_iterator<T, N>::diag_iterator(mm::square_matrix<T, N>& _M, int pos)
    index(0), M(_M), position(pos)
{
    assert(abs(pos) < N) // pos bounded between ]-N, N[ 
}

template<typename T, std::size_t N>
T& mm::diag_iterator<T, N>::operator*() const
{
    return (k > 0) ?
        M.data[(index + position) * Cols + index] :
        M.data[index * Cols + (index - position)];
}


template<typename T, std::size_t Rows, std::size_t Cols>
mm::const_vector_iterator<T, Rows, Cols>::const_vector_iterator(const mm::basic_matrix<T, Rows, Cols>& _M, std::size_t pos, bool dir)
    index(0), M(_M), position(pos), direction(dir)
{
    assert((dir && pos < Cols) || (!dir && pos < Rows))
}

template<typename T, std::size_t Rows, std::size_t Cols>
const T& mm::const_vector_iterator<T, Rows, Cols>::operator*() const
{
    return (direction) ?
            M.data[position * Cols + index] :
            M.data[index * Cols + position];
}

template<typename T, std::size_t Rows, std::size_t Cols>
const T& mm::const_vector_iterator<T, Rows, Cols>::operator[](std::size_t i) const
{
    return (direction) ?
            M.data[position * Cols + i] :
            M.data[i * Cols + position];
}

template<typename T, std::size_t N>
mm::const_diag_iterator<T, N>::const_diag_iterator(const mm::square_matrix<T, N>& _M, int pos)
    index(0), M(_M), position(pos)
{
    assert(abs(pos) < N) // pos bounded between ]-N, N[ 
}

template<typename T, std::size_t N>
T& mm::const_diag_iterator<T, N>::operator*() const
{
    return (k > 0) ?
        M.data[(index + position) * Cols + index] :
        M.data[index * Cols + (index - position)];
}


/* Iterators implementations */

template<typename T, std::size_t Rows, std::size_t Cols>
mm::vector_iterator<T, Rows, Cols>::vector_iterator(mm::basic_matrix<T, Rows, Cols>& _M, std::size_t pos, bool dir)
    index(0), M(_M), position(pos), direction(dir)
{
    assert((dir && pos < Cols) || (!dir && pos < Rows))
}

template<typename T, std::size_t Rows, std::size_t Cols>
T& mm::vector_iterator<T, Rows, Cols>::operator*() const
{
    return (direction) ?
            M.data[position * Cols + index] :
            M.data[index * Cols + position];
}

template<typename T, std::size_t Rows, std::size_t Cols>
T& mm::vector_iterator<T, Rows, Cols>::operator[](std::size_t i)
{
    return (direction) ?
            M.data[position * Cols + i] :
            M.data[i * Cols + position];
}

template<typename T, std::size_t N>
mm::diag_iterator<T, N>::diag_iterator(mm::square_matrix<T, N>& _M, int pos)
    index(0), M(_M), position(pos)
{
    assert(abs(pos) < N) // pos bounded between ]-N, N[ 
}

template<typename T, std::size_t N>
T& mm::diag_iterator<T, N>::operator*() const
{
    return (k > 0) ?
        M.data[(index + position) * Cols + index] :
        M.data[index * Cols + (index - position)];
}


template<typename T, std::size_t Rows, std::size_t Cols>
mm::const_vector_iterator<T, Rows, Cols>::const_vector_iterator(const mm::basic_matrix<T, Rows, Cols>& _M, std::size_t pos, bool dir)
    index(0), M(_M), position(pos), direction(dir)
{
    assert((dir && pos < Cols) || (!dir && pos < Rows))
}

template<typename T, std::size_t Rows, std::size_t Cols>
const T& mm::const_vector_iterator<T, Rows, Cols>::operator*() const
{
    return (direction) ?
            M.data[position * Cols + index] :
            M.data[index * Cols + position];
}

template<typename T, std::size_t Rows, std::size_t Cols>
const T& mm::const_vector_iterator<T, Rows, Cols>::operator[](std::size_t i) const
{
    return (direction) ?
            M.data[position * Cols + i] :
            M.data[i * Cols + position];
}

template<typename T, std::size_t N>
mm::const_diag_iterator<T, N>::const_diag_iterator(const mm::square_matrix<T, N>& _M, int pos)
    index(0), M(_M), position(pos)
{
    assert(abs(pos) < N) // pos bounded between ]-N, N[ 
}

template<typename T, std::size_t N>
T& mm::const_diag_iterator<T, N>::operator*() const
{
    return (k > 0) ?
        M.data[(index + position) * Cols + index] :
        M.data[index * Cols + (index - position)];
}


/* Iterators implementations */

template<typename T, std::size_t Rows, std::size_t Cols>
mm::vector_iterator<T, Rows, Cols>::vector_iterator(mm::basic_matrix<T, Rows, Cols>& _M, std::size_t pos, bool dir)
    index(0), M(_M), position(pos), direction(dir)
{
    assert((dir && pos < Cols) || (!dir && pos < Rows))
}

template<typename T, std::size_t Rows, std::size_t Cols>
T& mm::vector_iterator<T, Rows, Cols>::operator*() const
{
    return (direction) ?
            M.data[position * Cols + index] :
            M.data[index * Cols + position];
}

template<typename T, std::size_t Rows, std::size_t Cols>
T& mm::vector_iterator<T, Rows, Cols>::operator[](std::size_t i)
{
    return (direction) ?
            M.data[position * Cols + i] :
            M.data[i * Cols + position];
}

template<typename T, std::size_t N>
mm::diag_iterator<T, N>::diag_iterator(mm::square_matrix<T, N>& _M, int pos)
    index(0), M(_M), position(pos)
{
    assert(abs(pos) < N) // pos bounded between ]-N, N[ 
}

template<typename T, std::size_t N>
T& mm::diag_iterator<T, N>::operator*() const
{
    return (k > 0) ?
        M.data[(index + position) * Cols + index] :
        M.data[index * Cols + (index - position)];
}


template<typename T, std::size_t Rows, std::size_t Cols>
mm::const_vector_iterator<T, Rows, Cols>::const_vector_iterator(const mm::basic_matrix<T, Rows, Cols>& _M, std::size_t pos, bool dir)
    index(0), M(_M), position(pos), direction(dir)
{
    assert((dir && pos < Cols) || (!dir && pos < Rows))
}

template<typename T, std::size_t Rows, std::size_t Cols>
const T& mm::const_vector_iterator<T, Rows, Cols>::operator*() const
{
    return (direction) ?
            M.data[position * Cols + index] :
            M.data[index * Cols + position];
}

template<typename T, std::size_t Rows, std::size_t Cols>
const T& mm::const_vector_iterator<T, Rows, Cols>::operator[](std::size_t i) const
{
    return (direction) ?
            M.data[position * Cols + i] :
            M.data[i * Cols + position];
}

template<typename T, std::size_t N>
mm::const_diag_iterator<T, N>::const_diag_iterator(const mm::square_matrix<T, N>& _M, int pos)
    index(0), M(_M), position(pos)
{
    assert(abs(pos) < N) // pos bounded between ]-N, N[ 
}

template<typename T, std::size_t N>
T& mm::const_diag_iterator<T, N>::operator*() const
{
    return (k > 0) ?
        M.data[(index + position) * Cols + index] :
        M.data[index * Cols + (index - position)];
}


