#include "mm/mmmatrix.hpp"
#include "mm/view.hpp"

#include <iostream>
#include <complex>

// int main(int argc, char *argv[]) {
int main() {
	// std::cout << "MxN dimensional (int) matrices" << std::endl;

    mm::matrix<int, 2, 2> a { 1, 2, 3, 4 };

    // mm::matrix<int, 3, 2> a {{1, 2}, {3, 4}, {5, 6}};
    // mm::matrix<int, 3, 2> a {1, 2, 3, 4, 5, 6};

    // mm::matrix<int, 3, 2> b {4, 3, 9, 1, 2, 5};
    // mm::matrix<int, 2, 4> c {1, 2, 3, 4, 5, 6, 7, 8};

    // std::cout << "a = \n" << a;
    // std::cout << "b = \n" << b;
    // std::cout << "c = \n" << c;
    // std::cout << std::endl;
    
    // access elements
    // std::cout << "Access elements" << std::endl;
    // std::cout << "a.at(2,0) = " << a.at(2, 0) << std::endl;
    // std::cout << "a[2][0]   = " << a[2][0] << std::endl;;
    // std::cout << std::endl;
    
    // basic operations
    // std::cout << "Basic operations" << std::endl;
    // std::cout << "a + b = \n" << a + b;
    // std::cout << "a - b = \n" << a - b;
    // std::cout << "a * c = \n" << a * c;
    // std::cout << "a * 2 = \n" << a * 2;
    // std::cout << "2 * a = \n" << 2 * a;
    // std::cout << "a.td() = \n" << a.t(); // or a.trasposed();
    // std::cout << std::endl;

    std::cout << "a = \n" << a;

    std::cout << "Cloning a" << std::endl;
    decltype(a) e = mm::clone(a) | mm::alg::transpose();
    std::cout << "e = \n" << e;
    std::cout << std::endl;

    std::cout << "Mutating a" << std::endl;
    mm::mutate(a) | mm::alg::transpose();
    std::cout << "a = \n" << a;
    std::cout << std::endl;

    a | mm::alg::tr();
    std::cout << "a = \n" << a;

    // std::cout << "Converting clone object" << std::endl;
    // mm::matrix<int, 2, 2> g = e;
    // std::cout << std::endl;

    // std::cout << "Converting mutate object" << std::endl;
    // mm::matrix<int, 2, 2> h = f;
    // std::cout << std::endl;

	return 0;
}
