project(mlipp)
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
  add_compile_options(--std=c++1z -O0 -Wall)
else()
  add_compile_options(--std=c++1z -O3 -Wall)
endif()

include_directories(src/core)
include_directories(src/morton)

# add_executable(example_kdtree src/examples/example_kdtree.cpp)
# add_executable(example src/examples/example.cpp)
add_executable(example_zorder src/examples/example_zorder.cpp)
add_executable(example_mlipp src/examples/example_mlipp.cpp)
# add_executable(example_quad src/examples/example_quad.cpp)
# add_executable(example_without_fmcd src/examples/example_without_fmcd.cpp)
# add_executable(example_bulk_load src/examples/example_bulk_load.cpp)

# add_executable(point_1d_benchmark src/benchmarks/point_1d_benchmark.cpp)
# add_executable(point_benchmark src/benchmarks/point_benchmark.cpp)
# add_executable(range_benchmark src/benchmarks/range_benchmark.cpp)
# add_executable(range_benchmark_1d src/benchmarks/range_benchmark_1d.cpp)

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)


# Enable testing
enable_testing()
add_executable(
  point_test
  test/test_point.cc
)
target_link_libraries(
  point_test
  GTest::gtest_main
)

add_executable(
  mlipp_kd_test
  test/test_mlipp_kd.cc
)
target_link_libraries(
  mlipp_kd_test
  GTest::gtest_main
)

add_executable(
  test_mlipp_z
  test/test_mlipp_z.cc
)
target_link_libraries(
  test_mlipp_z
  GTest::gtest_main
)

add_executable(
  test_lippxy
  test/test_lippxy.cc
)
target_link_libraries(
  test_lippxy
  GTest::gtest_main
)

add_executable(
  test_multi_mlipp
  test/test_multi_mlipp.cc
)
target_link_libraries(
  test_multi_mlipp
  GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(mlipp_kd_test)
gtest_discover_tests(point_test)
gtest_discover_tests(test_mlipp_z)
gtest_discover_tests(test_lippxy)
gtest_discover_tests(test_multi_mlipp)
