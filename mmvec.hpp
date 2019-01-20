/* mmvec.hpp
 * Part of Mathematical library built (ab)using Modern C++ 17 abstractions.
 *
 * This library is not intended to be _performant_, it does not contain
 * hand written SMID / SSE / AVX optimizations. It is instead an example
 * of highly abstracted code, where Vectors can contain any data type.
 *
 * As a challenge, the vector data structure has been built on a container
 * of static capacity. But if a dynamic base container is needed, the code
 * should be easily modifiable to add further abstraction, by templating
 * the container, and by consequence the allocator.
 *
 * Naoki Pross <naopross@thearcway.org>
 * 2018 ~ 2019
 */
#include <iostream>

#include <cassert>
#include <cmath>

#include <array>
#include <algorithm>
#include <numeric>
#include <complex>
#include <initializer_list>

namespace mm {
    // generic implementation
    template<typename T, std::size_t d>
    struct basic_vec;

    // usable specializations
    template<typename T, std::size_t d>
    struct vec;
    template<typename T>
    struct vec3;
    template<typename T>
    struct vec2;
}

template<typename T, std::size_t d>
struct mm::basic_vec : public std::array<T, d> {
    using type = T;
    static constexpr std::size_t dimensions = d;

    // TODO: template away these
    static constexpr T null_element = static_cast<T>(0);
    static constexpr T unit_element = static_cast<T>(1);
    static constexpr T unit_additive_inverse_element = static_cast<T>(-1);

    basic_vec();
    basic_vec(const std::initializer_list<T> l);
    template<std::size_t n> basic_vec(const basic_vec<T, n>& other);

    T length() const;

    template<std::size_t n>
    basic_vec<T, d>& operator=(const mm::basic_vec<T, n>& other);

    template<std::size_t n>
    basic_vec<T, d>& operator+=(const mm::basic_vec<T, n>& other);

    template<std::size_t n>
    basic_vec<T, d>& operator-=(const mm::basic_vec<T, n>& other);


    basic_vec<T, d>& operator*=(const T& scalar);
};


// member functions for basic_vec

template<typename T, std::size_t d>
mm::basic_vec<T, d>::basic_vec() : std::array<T, d>() {
    this->fill(basic_vec<T, d>::null_element);
}

template<typename T, std::size_t d>
mm::basic_vec<T, d>::basic_vec(const std::initializer_list<T> l) {
    // construct with empty values
    basic_vec();

    // why can't this sh*t be a constexpr with static_assert???
    assert(l.size() <= d);
    std::copy(l.begin(), l.end(), this->begin());
}

template<typename T, std::size_t d>
template<std::size_t n>
mm::basic_vec<T, d>::basic_vec(const mm::basic_vec<T, n>& other) {
    // construct with empty values
    basic_vec();
    // uses operator=
    *this = other;
}

template<typename T, std::size_t d>
T mm::basic_vec<T, d>::length() const {
    return std::sqrt(std::accumulate(this->begin(), this->end(),
        basic_vec<T, d>::null_element,
        [](const T& init, const T& val) -> T {
            return init + val * val;
        }
    ));
}


// memeber operator overloads for basic_vec

template<typename T, std::size_t d>
template<std::size_t n>
mm::basic_vec<T, d>& mm::basic_vec<T, d>::operator=(const mm::basic_vec<T, n>& other) {
    static_assert(
        d >= n, 
        "cannot copy higher dimensional vector into a smaller one"
    );

    std::copy(other.begin(), other.end(), this->begin());

    return *this;
}

template<typename T, std::size_t d>
template<std::size_t n>
mm::basic_vec<T, d>& mm::basic_vec<T, d>::operator+=(const mm::basic_vec<T, n>& other) {
    *this = *this + other;
    return *this;
}

template<typename T, std::size_t d>
template<std::size_t n>
mm::basic_vec<T, d>& mm::basic_vec<T, d>::operator-=(const mm::basic_vec<T, n>& other) {
    *this = *this - other;
    return *this;
}

template<typename T, std::size_t d>
mm::basic_vec<T, d>& mm::basic_vec<T, d>::operator*=(const T& scalar) {
    *this = *this * scalar;
    return *this;
}


// operator overloads for basic_vec

template<typename T, std::size_t d>
mm::basic_vec<T, d> operator+(const mm::basic_vec<T, d>& rhs, const mm::basic_vec<T, d>& lhs) {
    mm::basic_vec<T, d> out;
    
    std::transform(rhs.begin(), rhs.end(), lhs.begin(), out.begin(),
        [](const T& r, const T& l) -> T {
            return r + l;
        }
    );

    return out;
}

template<typename T, std::size_t d>
mm::basic_vec<T, d> operator*(const mm::basic_vec<T, d>& rhs, const T& lhs) {
    return lhs * rhs;
}

template<typename T, std::size_t d>
mm::basic_vec<T, d> operator*(const T& rhs, const mm::basic_vec<T, d>& lhs) {
    mm::basic_vec<T, d> out;

    std::transform(lhs.begin(), lhs.end(), out.begin(), 
        [rhs](const T& t) -> T {
            return t * rhs;
    });

    return out;
}

template<typename T, std::size_t d>
mm::basic_vec<T, d> operator-(const mm::basic_vec<T, d>& rhs, const mm::basic_vec<T, d>& lhs) {
    return rhs + mm::basic_vec<T, d>::unit_additive_inverse_element * lhs;
}

template<typename T, std::size_t d>
T operator*(const mm::basic_vec<T, d>& rhs, const mm::basic_vec<T, d>& lhs) {
    return std::inner_product(rhs.begin(), rhs.end(), lhs.begin(), 0);
}

template<typename T, std::size_t d>
std::ostream& operator<<(std::ostream& os, const mm::basic_vec<T, d>& v) {
    os << "<";
    std::for_each(v.begin(), v.end() -1, [&](const T& el) {
        os << el << ", ";
    });
    os << v.back() << ">";

    return os;
}


// actual vectors to use in your code

template<typename T, std::size_t d>
class mm::vec: public mm::basic_vec<T, d> {
public:
    vec(std::initializer_list<T> l) : basic_vec<T, d>(l) {}

    template<std::size_t n>
    vec(const basic_vec<T, n>& other) : basic_vec<T, d>(other) {}
};


// three dimensional specialization with a static cross product
// TODO: specialize operator+ for spherical coordinates

template<typename T>
class mm::vec3 : public mm::basic_vec<T, 3> {
public:
    vec3() : basic_vec<T, 3>() {}
    vec3(std::initializer_list<T> l) : basic_vec<T, 3>(l) {}

    template<std::size_t n>
    vec3(const basic_vec<T, n>& other) : basic_vec<T, 3>(other) {}

    T& x() { return this->at(0); }
    T& y() { return this->at(1); }
    T& z() { return this->at(2); }

    const T& x() const { return this->at(0); }
    const T& y() const { return this->at(1); }
    const T& z() const { return this->at(2); }

    T zenith() const;
    T azimuth() const;
    vec3<T> spherical() const;

    static vec3<T> cross(const vec3<T>& rhs, const vec3<T>& lhs);
};

template<typename T>
T mm::vec3<T>::zenith() const {
    return std::acos(this->z() / this->length());
}

template<typename T>
T mm::vec3<T>::azimuth() const {
    return std::atan(this->y() / this->x());
}

template<typename T>
mm::vec3<T> mm::vec3<T>::spherical() const {
    return mm::vec3<T> {
        this->length(),
        this->zenith(),
        this->azimuth(),
    };
}

template<typename T>
mm::vec3<T> mm::vec3<T>::cross(const vec3<T>& rhs, const vec3<T>& lhs) {
    mm::vec3<T> res;

    res.x() = (rhs.y() * lhs.z()) - (rhs.z() * lhs.y());
    res.y() = (rhs.z() * lhs.x()) - (rhs.x() * lhs.z());
    res.z() = (rhs.x() * lhs.y()) - (rhs.y() * lhs.x());

    return res;
}


// two dimensional specialization with a polar conversion
// TODO: specialize operator+ for polar coordinates

template<typename T>
class mm::vec2: public mm::basic_vec<T, 2> {
public:
    vec2() : basic_vec<T, 2>() {}
    vec2(std::initializer_list<T> l) : basic_vec<T, 2>(l) {}

    template<std::size_t n>
    vec2(const basic_vec<T, n>& other) : basic_vec<T, 2>(other) {}

    T& x() { return this->at(0); }
    T& y() { return this->at(1); }

    const T& x() const { return this->at(0); }
    const T& y() const { return this->at(1); }

    T angle() const;
    vec2<T> polar() const;

    static vec3<T> cross(const vec2<T>& rhs, const vec2<T>& lhs);
};

template<typename T>
T mm::vec2<T>::angle() const {
    return std::atan(this->y() / this->x());
}

template<typename T>
mm::vec2<T> mm::vec2<T>::polar() const {
    return mm::vec2 {
        this->length(),
        this->angle()
    };
}

template<typename T>
mm::vec3<T> mm::vec2<T>::cross(const mm::vec2<T>& rhs, const mm::vec2<T>& lhs) {
    return mm::vec3<T>::cross(mm::vec3<T>(rhs), mm::vec3<T>(lhs));
}
