cmake_minimum_required (VERSION 2.6)
project (GeoLib)

set(EXECUTABLE_OUTPUT_PATH build)

add_executable(geloib_tests tests/geodistance_test.cpp)

target_compile_options(geloib_tests PUBLIC -std=c++11 -lboost_serialization -larmadillo -lmlpack -DBOOST_TEST_DYN_LINK)
