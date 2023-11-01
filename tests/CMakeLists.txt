enable_testing()

add_executable(main-test main-test.cpp)
add_test(NAME main-test COMMAND main-test)

target_link_libraries(main-test PRIVATE butterworth -lsndfile)

add_executable(test-butter test-butter.cpp)
add_test(NAME test-butter COMMAND test-butter)

