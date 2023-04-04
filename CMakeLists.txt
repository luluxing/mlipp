project(lipp)
cmake_minimum_required(VERSION 2.9)

set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Comment out code used for debugging purposes so it is not concerned by the coverage
if(CMAKE_BUILD_TYPE MATCHES Debug)
  message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
  add_definitions(-DDEBUG)
endif()

add_compile_options(--std=c++11 -O3 -Wall)

include_directories(src/core)
include_directories(src/morton)

add_executable(example_kdtree src/examples/example_kdtree.cpp)
add_executable(example src/examples/example.cpp)
add_executable(example_zorder src/examples/example_zorder.cpp)
add_executable(example_mlipp src/examples/example_mlipp.cpp)
add_executable(example_quad src/examples/example_quad.cpp)
add_executable(example_without_fmcd src/examples/example_without_fmcd.cpp)
add_executable(example_bulk_load src/examples/example_bulk_load.cpp)
