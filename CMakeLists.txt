cmake_minimum_required(VERSION 3.0.2)
project(lab1)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

find_package (Eigen3 3.4.0 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

add_executable(standar_multi src/standard_matrix_multiplication_gtod.cpp)
add_executable(eigen_multi src/eigen_matrix_multiplication_gtod.cpp)
