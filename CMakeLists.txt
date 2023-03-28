project(lipp)
cmake_minimum_required(VERSION 2.9)

add_compile_options(--std=c++11 -O3 -Wall)

include_directories(src/core)

add_executable(example_kd src/examples/kdtree.cpp)
add_executable(example src/examples/example.cpp)
add_executable(example_2d src/examples/example_2d.cpp)
add_executable(example_2darray src/examples/example_2darray.cpp)
add_executable(example_without_fmcd src/examples/example_without_fmcd.cpp)
add_executable(example_bulk_load src/examples/example_bulk_load.cpp)
