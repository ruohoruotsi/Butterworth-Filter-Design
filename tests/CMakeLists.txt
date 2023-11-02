enable_testing()

if(NOT ARM64)
    add_executable(main-test main-test.cpp)
    add_test(NAME main-test COMMAND main-test)
    target_link_libraries(main-test PRIVATE butterworth -lsndfile)
endif()

add_executable(test-butter test-butter.cpp)
add_test(NAME test-butter COMMAND test-butter)

add_executable(test-conv test-conv.cpp)
add_test(NAME test-conv COMMAND test-conv)
target_link_libraries(test-conv PRIVATE butterworth)
