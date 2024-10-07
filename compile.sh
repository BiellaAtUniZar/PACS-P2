#!/bin/bash

# Define the file paths
STANDARD_FILE="standard_matrix_multiplication.cpp"
EIGEN_FILE="eigen_matrix_multiplication.cpp"
EIGEN_INCLUDE_DIR="./vendor/eigen-3.4.0/"

# Define output files
OUT_STANDARD_NO_OPT="standard_matrix_multiplication_no-o"
OUT_STANDARD_O2="standard_matrix_multiplication_o2"
OUT_STANDARD_O3="standard_matrix_multiplication_o3"

OUT_EIGEN_NO_OPT="eigen_matrix_multiplication_no-o"
OUT_EIGEN_O2="eigen_matrix_multiplication_o2"
OUT_EIGEN_O3="eigen_matrix_multiplication_o3"

# Compile without optimization
g++ -o $OUT_STANDARD_NO_OPT $STANDARD_FILE
g++ -o $OUT_EIGEN_NO_OPT -I $EIGEN_INCLUDE_DIR $EIGEN_FILE

# Compile with O2 optimization
g++ -O2 -o $OUT_STANDARD_O2 $STANDARD_FILE
g++ -O2 -o $OUT_EIGEN_O2 -I $EIGEN_INCLUDE_DIR $EIGEN_FILE

# Compile with O3 optimization
g++ -O3 -o $OUT_STANDARD_O3 $STANDARD_FILE
g++ -O3 -o $OUT_EIGEN_O3 -I $EIGEN_INCLUDE_DIR $EIGEN_FILE

echo "Compilation done!"
