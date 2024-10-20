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


int main(int argc, char *argv[]) {

	int SIZE;
	bool AUTO_MODE = argc != 3;
	int iteration;
	if (AUTO_MODE) {
		iteration = 10;
		SIZE = 10;
	} else {
		SIZE = stoi(argv[1]); // Parses parameter, throws exception if not int
		iteration = stoi(argv[2]);
	}
	
	float mean_setup_time=0;
	std::vector<float> setup_times;
	
	for(int i=0;i<iteration;i++){

	timeval setup_start, setup_end;
	gettimeofday(&setup_start, NULL);
	// Fills matrices
	Mat m1 = fill_matrix(SIZE);
	Mat m2 = fill_matrix(SIZE);
	gettimeofday(&setup_end, NULL);
	// setup time
	auto ssec_diff = setup_end.tv_sec*1.+setup_end.tv_usec*(1e-6)-(setup_start.tv_sec*1.+setup_start.tv_usec*(1e-6));
	ssec_diff*=1000.;
	mean_setup_time+=ssec_diff;
	setup_times.push_back(ssec_diff);
	}
	
	mean_setup_time /= iteration * 1.;
	double stand_dev_setup = 0.;
	for (size_t i = 0; i < setup_times.size(); i++) {
		stand_dev_setup += (setup_times[i] - mean_setup_time) * (setup_times[i] - mean_setup_time) / setup_times.size() * 1.0;
	}
	stand_dev_setup = sqrt(stand_dev_setup);
	std::cout << "Setup= " << mean_setup_time << "ms ; Standard Deviation= " << stand_dev_setup
						<< "ms" << std::endl;

	// matrix multiplication
	// measure performance time for mat mul
	float mean_exec_time=0;
		// Fills matrices
	Mat m1 = fill_matrix(SIZE);
	Mat m2 = fill_matrix(SIZE);
	std::vector<float> exec_times;
	for(int i=0;i<iteration;i++){
		timeval t1, t2;
		gettimeofday(&t1, NULL);
		Mat result = mul_matrix(m1, m2);
		gettimeofday(&t2, NULL);
		auto sec_diff = t2.tv_sec*1.+t2.tv_usec*(1e-6)-(t1.tv_sec*1.+t1.tv_usec*(1e-6));
		sec_diff*=1000.;
		mean_exec_time+=sec_diff;
		exec_times.push_back(sec_diff);
	}


	mean_exec_time /= iteration * 1.;
	double stand_dev = 0.;
	for (size_t i = 0; i < exec_times.size(); i++) {
		stand_dev +=
				(exec_times[i] - mean_exec_time) * (exec_times[i] - mean_exec_time) / exec_times.size() * 1.0;
	}
	stand_dev = sqrt(stand_dev);
	std::cout << "Exec= " << mean_exec_time << "ms ; Standard Deviation= " << stand_dev
						<< "ms" << std::endl;

	return 0;
}
