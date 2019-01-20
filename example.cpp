#include "mmvec.hpp"

#include <iostream>
#include <complex>

int main(int argc, char *argv[]) {
	// N dimensional vectors
	std::cout << "N dimensional (int) vectors" << std::endl;
	mm::vec<int, 5> u {3, 2, 1, 0, 1};
	mm::vec<int, 5> v {1, 2, 3, 4, 5};

	std::cout << "u = " << u << std::endl;
	std::cout << "v = " << v << std::endl;
	std::cout << std::endl;

	// basic operations
	std::cout << "u + v = " << u + v << std::endl;
	std::cout << "u - v = " << u - v << std::endl;
	std::cout << "2 * v = " << 2 * v << std::endl;
	std::cout << "v * 2 = " << v * 2 << std::endl;
	std::cout << "u * v = " << u * v << std::endl;
	std::cout << std::endl;

	// three dimensional vectors
	std::cout << "three dimensional (double) vectors" << std::endl;

	mm::vec3<double> a {1, 2, 3};
	mm::vec3<double> b {3, 2, 1};

	std::cout << "a = " << a << std::endl;
	std::cout << "b = " << b << std::endl;
	std::cout << std::endl;

	std::cout << "a x b        = " << mm::vec3<double>::cross(a, b) << std::endl;
	std::cout << "zenith(a)    = " << a.zenith() << std::endl;
	std::cout << "azimuth(a)   = " << a.azimuth() << std::endl;
	std::cout << "spherical(a) = " << a.spherical() << std::endl;
	std::cout << std::endl;

	// two dimensional vector
	std::cout << "two dimensional (complex) vectors" << std::endl;

	mm::vec2<std::complex<float>> j {{1, 2}, {3, -1}};
	mm::vec2<std::complex<float>> k { 5, {-2, 1}};

	std::cout << "j = " << j << std::endl;
	std::cout << "k = " << k << std::endl;
	std::cout << std::endl;

	std::cout << "j x k    = " << mm::vec2<std::complex<float>>::cross(j, k) << std::endl;
	std::cout << "angle(j) = " << j.angle() << std::endl;
	std::cout << "polar(j) = " << j.polar() << std::endl;
	std::cout << std::endl;

	return 0;
}
