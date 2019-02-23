add_executable(vec_example
    vec_example.cpp
)

target_compile_options(vec_example
    PRIVATE
        $<$<CXX_COMPILER_ID:GNU>:
            -pedantic -Wall -Wextra -Wcast-qual -Wcast-align -Wpointer-arith
            -Winit-self -Wshadow -Wswitch-enum -Wredundant-decls -Wfloat-equal
            -Wundef -Wvla -Wconversion -Wstrict-aliasing
        >
        $<$<CXX_COMPILER_ID:MSVC>:/W4>
)

target_compile_features(vec_example
    PRIVATE
        cxx_std_17
)

target_include_directories(vec_example
    PRIVATE
        ${MM_INCLUDE_DIRS}
)

target_link_libraries(vec_example
    PRIVATE
        MM::mm
)
