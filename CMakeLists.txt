cmake_minimum_required(VERSION 3.0)
project(crg)


list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")


file(STRINGS config/version.txt GALOIS_VERSION)
string(REGEX REPLACE "[ \t\n]" "" GALOIS_VERSION ${GALOIS_VERSION})
string(REGEX REPLACE "([0-9]+)\\.([0-9]+)\\.([0-9]+)" "\\1" GALOIS_VERSION_MAJOR ${GALOIS_VERSION})
string(REGEX REPLACE "([0-9]+)\\.([0-9]+)\\.([0-9]+)" "\\2" GALOIS_VERSION_MINOR ${GALOIS_VERSION})
string(REGEX REPLACE "([0-9]+)\\.([0-9]+)\\.([0-9]+)" "\\3" GALOIS_VERSION_PATCH ${GALOIS_VERSION})
set(GALOIS_COPYRIGHT_YEAR "2018") # Also in COPYRIGHT
set(CMAKE_BUILD_TYPE "Release")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_COMPILER /usr/bin/g++-11)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
SET (CMAKE_CXX_FLAGS "-g -O3 -fopenmp -shared-libgcc -fcoroutines")


find_package(Threads REQUIRED)

find_package(Boost 1.58.0 REQUIRED COMPONENTS serialization iostreams)


add_custom_target(lib)
add_custom_target(apps)

add_subdirectory(libcoro)
add_subdirectory(app)
