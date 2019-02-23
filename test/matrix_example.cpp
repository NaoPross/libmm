#include "mmmatrix.hpp"

#include <iostream>
#include <complex>

int main(int argc, char *argv[]) {
	std::cout << "MxN dimensional (int) matrices" << std::endl;
    mm::matrix<int, 3, 2> m {{1, 2}, {3, 4}, {5, 6}};

    std::cout << m;

	return 0;
}
