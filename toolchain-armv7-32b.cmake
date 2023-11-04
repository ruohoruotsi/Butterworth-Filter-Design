set(CMAKE_SYSTEM_NAME Linux)                          # Operating system type
set(CMAKE_SYSTEM_PROCESSOR arm)                       # CPU architecture

set(CMAKE_C_COMPILER arm-linux-musleabihf-gcc)        # gcc compiler program name
set(CMAKE_CXX_COMPILER arm-linux-musleabihf-g++)      # g++ compiler program name

set(CMAKE_FIND_ROOT_PATH /opt/cross/)          # Root directory of g++
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM BOTH)    # Tell CMake to look for programs like g++ from under ROOT_PATH, but it can also look elsewhere
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)    # Tell CMake to look for library files from under ROOT_PATH
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)    # Tell CMake to look for header files from under ROOT_PATH

set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_C_COMPILER_TARGET arm-linux-musleabihf)    # Target programs that use the musl libc library, if they are standard gnu libraries, can comment out this line

# set(CMAKE_C_FLAGS "-static-libgcc -static-libstdc++")

# Set compilation options
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=armv7-a -mfpu=vfpv3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=armv7-a -mfpu=vfpv3")    # fpu版本

# Set link options
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath-link,/opt/cross/arm-linux-musleabihf/lib")
