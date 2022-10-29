cmake_minimum_required(VERSION 3.21)

set(CMAKE_OSX_ARCHITECTURES x86_64)
project(RCI C)
set(CMAKE_C_STANDARD 11)
SET(CMAKE_BUILD_TYPE "Debug")

find_package(Matlab)
add_compile_definitions(MATLAB_MEX_FILE)

# Configure MKL paths
set(ENV{MKL_LIB_PATH} /opt/intel/oneapi/mkl/2021.2.0/lib)
set(ENV{MKL_OMP_PATH} /opt/intel/oneapi/compiler/2021.2.0/mac/compiler/lib)

# Include project directory
include_directories(.)
include_directories(/opt/intel/oneapi/mkl/2021.2.0/include)

matlab_add_mex(NAME rci SRC ./rci_cg.c ./matrix.h ./tmwtypes.h ./mex.h)

target_link_libraries(rci $ENV{MKL_LIB_PATH}/libmkl_core.a)
target_link_libraries(rci $ENV{MKL_LIB_PATH}/libmkl_intel_lp64.a)
target_link_libraries(rci $ENV{MKL_LIB_PATH}/libmkl_sequential.a)