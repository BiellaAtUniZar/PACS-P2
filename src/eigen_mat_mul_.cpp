#include <sys/time.h>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <random>
#include <cmath>
#include <vector>

using namespace std;

typedef Eigen::MatrixXd Mat;
/**
 * Fills the square matrix of dimension SIZExSIZE of random values between -10
 * and 10 included
 * @param matrix one dimensional pointer to the matrix
 * @param SIZE size of the matrix
 */
Mat fill_matrix(const size_t SIZE) {
	Eigen::MatrixXd mat = Eigen::MatrixXd::Random(SIZE, SIZE);
	return mat;
}

/**
 * Performs multiplication and returns address of resulting matrix
 * @param matrix1
 * @param matrix2
 * @param SIZE
 */
Mat mul_matrix(Mat mat1, Mat mat2) { return mat1 * mat2; }

void print_matrix(Mat mat1) { std::cout << mat1 << std::endl; }

/**
 * Function to get time difference in milliseconds
 */
double get_time_diff(const timeval& start, const timeval& end) {
		return (end.tv_sec - start.tv_sec) * 1000.0 + (end.tv_usec - start.tv_usec) / 1000.0;
}

int main(int argc, char *argv[]) {
	//Begin instrumentation
 
  	//End instrumentation


	int SIZE;
	bool AUTO_MODE = argc != 2;
	if (AUTO_MODE) {
		SIZE = 10;
	} else {
		SIZE = stoi(argv[1]); // Parses parameter, throws exception if not int
	}
	// Fills matrices
	Mat m1 = fill_matrix(SIZE);
	Mat m2 = fill_matrix(SIZE);
	///matrix multiplication
    Mat result = mul_matrix(m1, m2);


	return 0;
}
