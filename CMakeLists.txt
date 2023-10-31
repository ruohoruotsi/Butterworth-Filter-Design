cmake_minimum_required(VERSION 3.5)
project(butter_worth_filter_design)

include_directories(${CMAKE_SOURCE_DIR})

# lib_comfort库
add_library(butterworth-filter-design SHARED 
    Butterworth.cpp
    Biquad.cpp
  )

# 添加测试
enable_testing()
add_subdirectory(tests)
