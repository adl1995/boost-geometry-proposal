cmake_minimum_required (VERSION 2.6)
project (GeoLib)

SET(CMAKE_CXX_COMPILER /usr/bin/g++-7)

set(EXECUTABLE_OUTPUT_PATH bin)
add_executable(geloib_tests tests/geodistance_test.cpp)

target_link_libraries(geloib_tests -std=c++11 -lboost_serialization -lboost_system -lboost_timer -lboost_chrono)

add_custom_target(tests COMMAND ./bin/geloib_tests)

add_custom_target(docs COMMAND doxygen ../Doxyfile)
