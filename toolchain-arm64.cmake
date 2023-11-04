set(CMAKE_SYSTEM_NAME Linux)                          # Operating system type
set(CMAKE_SYSTEM_PROCESSOR aarch64)                     # CPU architecture
set(CMAKE_C_COMPILER_TARGET aarch64-linux-gnu)

set(CMAKE_C_COMPILER aarch64-linux-gnu-gcc)        # gcc compiler program name
set(CMAKE_CXX_COMPILER aarch64-linux-gnu-g++)      # g++ compiler program name
set(CMAKE_Fortran_COMPILER aarch64-linux-gnu-gfortran)      # g++ compiler program name  

set(CMAKE_FIND_ROOT_PATH /opt/cross-gcc-linaro-6.3.1-arm64/)          # Root directory of g++
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM BOTH)    # Tell CMake to look for programs like g++ from under ROOT_PATH, but it can also look elsewhere
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)    # Tell CMake to look for library files from under ROOT_PATH
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)    # Tell CMake to look for header files from under ROOT_PATH
# set(CMAKE_PREFIX_PATH /opt/cross-gcc-linaro-6.3.1-arm64/arm-linux-gnueabihf)

# Whether to statically link libc and libstdc++
# set(CMAKE_C_FLAGS "-static-libgcc -static-libstdc++")

# Set compilation options
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=armv7-a -mfpu=vfpv3")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=armv7-a -mfpu=vfpv3")    # fpu版本

# Set link options
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath-link,/opt/cross-gcc-linaro-6.3.1-arm64/aarch64-linux-gnu/lib")
