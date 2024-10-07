cmake_minimum_required(VERSION 3.0.2)
project(lab1)

find_package (Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
add_executable(standar_multi src/standard_matrix_multiplication_gtod.cpp)
add_executable(eigen_multi src/eigen_matrix_multiplication_gtod.cpp)
