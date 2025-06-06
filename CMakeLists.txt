## CMakeLists.txt
cmake_minimum_required(VERSION 3.18)
project("butterworth_fdatool" LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

# enable testing with CTest
include(CTest)

find_package(Catch2 3 REQUIRED)
include(Catch)

# These tests can use the Catch2-provided main
add_executable(tests tests.cpp Butterworth.cpp Biquad.cpp)

find_library(sndfile libsndfile)
target_link_libraries(tests PRIVATE Catch2::Catch2WithMain sndfile)
catch_discover_tests(tests)