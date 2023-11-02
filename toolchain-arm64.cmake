set(CMAKE_SYSTEM_NAME Linux)                          # 操作系统类型
set(CMAKE_SYSTEM_PROCESSOR aarch64)                     # CPU 架构
set(CMAKE_C_COMPILER_TARGET aarch64-linux-gnu)

set(CMAKE_C_COMPILER aarch64-linux-gnu-gcc)        # gcc 编译器程序名
set(CMAKE_CXX_COMPILER aarch64-linux-gnu-g++)      # g++ 编译器程序名  
set(CMAKE_Fortran_COMPILER aarch64-linux-gnu-gfortran)      # g++ 编译器程序名  

set(CMAKE_FIND_ROOT_PATH /opt/cross-gcc-linaro-6.3.1-arm64/)          # g++ 所在根目录
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM BOTH)    # 告诉CMake要从 ROOT_PATH 下查找g++等程序，但也可以从其他地方查找
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)    # 告诉CMake要从 ROOT_PATH 下查找library库文件
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)    # 告诉CMake要从 ROOT_PATH 下查找头文件
# set(CMAKE_PREFIX_PATH /opt/cross-gcc-linaro-6.3.1-arm64/arm-linux-gnueabihf)

# 是否静态链接 libc 和 libstdc++
# set(CMAKE_C_FLAGS "-static-libgcc -static-libstdc++")

# 设置编译选项
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=armv7-a -mfpu=vfpv3")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=armv7-a -mfpu=vfpv3")    # fpu版本

# 设置链接选项
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath-link,/opt/cross-gcc-linaro-6.3.1-arm64/aarch64-linux-gnu/lib")
