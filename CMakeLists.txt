# 交叉编译为 arm 架构程序。注意: 本块必须放在最前面，否则不生效
if(ARM64)
  message(STATUS "Using ARM64 Tool Chain!")
  set(CMAKE_TOOLCHAIN_FILE ${CMAKE_SOURCE_DIR}/toolchain-arm64.cmake)
  # 总是 Release 版本
  set(CMAKE_BUILD_TYPE Release)
  message(STATUS "Compiler:${CMAKE_CXX_COMPILER}")
else()
  message(STATUS "Using X64 Tool Chain!")
  # 指定编译器选项
  # 为了 gdb 调试，不优化、带上debug符号
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g -rdynamic")
endif()

cmake_minimum_required(VERSION 3.5)
project(butter_worth_filter_design)

set(CMAKE_TOOLCHAIN_FILE ${CMAKE_SOURCE_DIR}/toolchain-arm64.cmake)


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

# 添加测试
enable_testing()
add_subdirectory(tests)
