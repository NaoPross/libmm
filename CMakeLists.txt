cmake_minimum_required(VERSION 3.10)
project(
    MM
    VERSION 0.1
    DESCRIPTION "MiniMath, a mathematical library that (ab)uses abstraction"
    LANGUAGES CXX
)

list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_SOURCE_DIR}/cmake)

############################
# build libmm
add_library(mm STATIC
    ${CMAKE_CURRENT_SOURCE_DIR}/mmvec.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mmmatrix.cpp
)

add_library(MM::mm ALIAS mm)

# set up headers
target_include_directories(mm
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    PRIVATE
        ${CMAKE_CURRENT_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_compile_options(mm
    PRIVATE
        $<$<CXX_COMPILER_ID:GNU>:
            -pedantic -Wall -Wextra -Wcast-qual -Wcast-align -Wpointer-arith
            -Winit-self -Wshadow -Wswitch-enum -Wredundant-decls -Wfloat-equal
            -Wundef -Wvla -Wconversion -Wstrict-aliasing
        >
        $<$<CXX_COMPILER_ID:MSVC>:/W4>
)

target_compile_definitions(mm
    PRIVATE
        $<$<CONFIG:Debug>:DEBUG>
)

# set compiler features
# For a complete list see:
# https://cmake.org/cmake/help/latest/prop_gbl/CMAKE_CXX_KNOWN_FEATURES.html
target_compile_features(mm
    PRIVATE
        cxx_std_17
)


############################
# installation
include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/MM)

install(TARGETS mm
    EXPORT mm-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

set_target_properties(mm PROPERTIES EXPORT_NAME MM)


# install public headers
install(
    FILES 
        ${CMAKE_CURRENT_SOURCE_DIR}/include/mm/mmvec.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/mm/mm.hpp
    DESTINATION
        ${CMAKE_INSTALL_INCLUDEDIR}/mm
)

# export targets
install(EXPORT mm-targets
    FILE
        MMTargets.cmake
    NAMESPACE
        MM::
    DESTINATION
        ${INSTALL_CONFIGDIR}
)

# create ConfigVersion.cmake
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/MMConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(
    ${CMAKE_CURRENT_LIST_DIR}/cmake/MMConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/MMConfig.cmake
    INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)

# install config, configversion
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/MMConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/MMConfigVersion.cmake
    DESTINATION ${INSTALL_CONFIGDIR}
)

############################
# export from the build tree
export(EXPORT mm-targets
    FILE ${CMAKE_CURRENT_BINARY_DIR}/MMTargets.cmake
    NAMESPACE MM::
)

export(PACKAGE MM)

############################
# tests
add_subdirectory(test)
