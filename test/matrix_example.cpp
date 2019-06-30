#include "mm/mmmatrix.hpp"

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
    std::cout << std::endl; 
    
    // access elements
    std::cout << "Access elements" << std::endl;
    std::cout << "a.at(2,0) = " << a.at(2, 0) << std::endl;
    std::cout << "a[2][0]   = " << a[2][0] << std::endl;;
    std::cout << std::endl;
    
    /*
    // basic operations
    std::cout << "Basic operations" << std::endl;
    std::cout << "a + b = \n" << a + b;
    std::cout << "a - b = \n" << a - b;
    std::cout << "a * c = \n" << a * c;
    std::cout << "a * 2 = \n" << a * 2;
    std::cout << "2 * a = \n" << 2 * a;
    std::cout << "a.td() = \n" << a.td(); // or a.trasposed();
    std::cout << std::endl;*/

    // special matrices
    /*mm::square_matrix<std::complex<int>, 2> f {{{2, 3}, {1, 4}}, {{6, 1}, {-3, 4}}};

    std::cout << "Square matrix" << std::endl;
    std::cout << "f = \n" << f;

    std::cout << "tr(f) = " << f.tr();  //or f.trace()  << std::endl;

    mm::t_square_matrix<std::complex<int>, 2>& ft = f.t();
    std::cout << "after in place transpose f.t(), f = \n" << ft; 
    std::cout << std::endl;

    auto identity = mm::square_matrix<int, 3>::identity();

    std::cout << "Identity matrix" << std::endl;
    std::cout << "I = \n" << identity;
    std::cout << std::endl; */


	return 0;
}
