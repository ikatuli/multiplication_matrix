# cmake -DCMAKE_CXX_FLAGS:STRING="${CFLAGS}" ./

cmake_minimum_required(VERSION 3.12.3)
project(multiplication_matrix LANGUAGES CXX)
find_package(MPI REQUIRED) #find_package(MPI)
include_directories(${PROJECT_SOURCE_DIR})
#include_directories(SYSTEM ${MPI_INCLUDE_PATH})
add_executable(multiplication_matrix multiplication_matrix.cpp)
target_link_libraries(multiplication_matrix ${CMAKE_DL_LIBS} MPI::MPI_CXX)
set(CMAKE_BUILD_TYPE Debug)
