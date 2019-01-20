CPPC    := c++
CPPARGS := -Wall -Werror -I. -std=c++17 -pedantic

all: build/example

build/%: %.cpp mmvec.hpp
	mkdir -p build
	$(CPPC) $(CPPARGS) $< -o $@
