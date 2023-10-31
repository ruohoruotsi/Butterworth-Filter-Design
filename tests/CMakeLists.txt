enable_testing()

add_executable(test-butter test-butter.cpp)
add_test(NAME test-butter COMMAND test-butter)

