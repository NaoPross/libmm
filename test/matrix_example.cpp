#include "mm/mmmatrix.hpp"

#include <iostream>
#include <complex>

int main(int argc, char *argv[]) {
	std::cout << "MxN dimensional (int) matrices" << std::endl;
    mm::matrix<int, 3, 2> a {{1, 2}, {3, 4}, {5, 6}};
    mm::matrix<int, 3, 2> b {{4, 3}, {9, 1}, {2, 5}};
    mm::matrix<int, 2, 4> c {{1, 2, 3, 4}, {5, 6, 7, 8}};
    auto ct = c.t();

    std::cout << "a = \n" << a;
    std::cout << "b = \n" << b;
    std::cout << "c = \n" << c;
    std::cout << "c^t = \n" << ct;
    std::cout << std::endl; 
    
    // access elements
    std::cout << "Access elements" << std::endl;
    std::cout << "a.at(2,0) = " << a.at(2, 0) << std::endl;
    std::cout << "a[2][0]   = " << a[2][0] << std::endl;;
    std::cout << std::endl;
    
    // basic operations
    std::cout << "Basic operations" << std::endl;
    std::cout << "a + b = \n" << a + b;
    std::cout << "a - b = \n" << a - b;
    std::cout << "a * c = \n" << a * c;
    std::cout << "a * 2 = \n" << a * 2;
    std::cout << "2 * a = \n" << 2 * a;
    std::cout << "a.td() = \n" << a.t(); // or a.trasposed();
    std::cout << std::endl;

	return 0;
}
