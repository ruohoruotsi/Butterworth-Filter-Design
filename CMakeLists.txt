cmake_minimum_required(VERSION 3.5)
project(butter_worth_filter_design)

include_directories(${CMAKE_SOURCE_DIR})

# lib_comfort库
add_library(butterworth SHARED 
    Butterworth.cpp
    Biquad.cpp
    conv.cpp
    sos2tf.cpp
  )
# 添加编译选项，按位置无关代码编译
target_compile_options(butterworth PRIVATE -fPIC)

# 编译测试主程序 main
add_executable(test-main main.cpp)
target_link_libraries(test-main PRIVATE butterworth -lsndfile)

# 添加测试
enable_testing()
add_subdirectory(tests)
