#include "mmmatrix.hpp"

#include <iostream>
#include <complex>

int main(int argc, char *argv[]) {
	std::cout << "MxN dimensional (int) matrices" << std::endl;
    mm::matrix<int, 3, 2> a {{1, 2}, {3, 4}, {5, 6}};
    mm::matrix<int, 3, 2> b {{4, 3}, {9, 1}, {2, 5}};
    mm::matrix<int, 2, 4> c {{1, 2, 3, 4}, {5, 6, 7, 8}};

    std::cout << "a = \n" << a;
    std::cout << "b = \n" << b;
    std::cout << "c = \n" << c;

    // access elements
    std::cout << "a.at(2,0) = " << a.at(2, 0) << std::endl;
    std::cout << "a[2][0]   = " << a[2][0] << std::endl;;
    
    // basic operations
    std::cout << "a + b = \n" << a + b;
    std::cout << "a - b = \n" << a - b;
    std::cout << "a * c = \n" << a * c;
    std::cout << "a * 2 = \n" << a * 2;
    std::cout << "2 * a = \n" << 2 * a;
    std::cout << "tr(a) = \n" << a.trd();

	return 0;
}
