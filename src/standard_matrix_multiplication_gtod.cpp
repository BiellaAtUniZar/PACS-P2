#include <assert.h>
#include <chrono>
#include <iostream>
#include <random>
#include <vector>
#include <sys/time.h>

using namespace std;
using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;

typedef std::vector<std::vector<double>> Mat;
typedef std::vector<double> Row;
/**
 * Fills the square matrix of dimension SIZExSIZE of random values between -10
 * and 10 included
 * @param matrix one dimensional pointer to the matrix
 */
void fill_matrix(Mat &matrix)
{
	random_device rd; // Setup of number generator
	mt19937 gen(rd());
	uniform_real_distribution<> double_dist(-10.0, 10.0);
	int k = 0;
	for (size_t i = 0; i < matrix.size(); i++)
	{
		for (size_t j = 0; j < matrix[0].size(); j++)
		{
			// matrix[i][j] = double_dist(gen);  // Filling the array with random
			// doubles
			k++;
			matrix[i][j] = k;
		}
	}
}

/**
 * Performs multiplication and returns address of resulting matrix
 * @param matrix1
 * @param matrix2
 */
void mul_matrix(const Mat &matrix1, const Mat &matrix2, Mat &result)
{
	// Assuming NxN matrices
	assert(matrix1.size() == matrix2.size());
	assert(matrix1[0].size() == matrix2[0].size());
	assert(matrix1[0].size() == matrix1.size());
	assert(matrix2[0].size() == matrix2.size());
	for (size_t i = 0; i < matrix1.size(); i++)
	{
		for (size_t j = 0; j < matrix1.size(); j++)
		{
			for (size_t k = 0; k < matrix1.size(); k++)
			{
				result[i][j] += matrix1[i][k] * matrix2[k][j];
			}
		}
	}
}

void print_matrix(const Mat &matrix)
{
	for (size_t i = 0; i < matrix.size(); i++)
	{
		for (size_t j = 0; j < matrix[0].size(); j++)
		{
			std::cout << matrix[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

int main(int argc, char *argv[])
{

	// Begin instrumentation

	// End instrumentation

	int SIZE;
	bool AUTO_MODE = argc != 3;
	int iteration;
	if (AUTO_MODE)
	{
		iteration = 10;
		SIZE = 10;
	}
	else
	{
		SIZE = stoi(argv[1]); // Parses parameter, throws exception if not int
		iteration = stoi(argv[2]);
	}

	float mean_setup_time = 0;
	std::vector<float> setup_times;
	for (int i = 0; i < iteration; i++)
	{

		timeval setup_start, setup_end;
		gettimeofday(&setup_start, NULL);
		// Fills matrices
		Mat m1(SIZE, Row(SIZE, 0.));
		Mat m2(SIZE, Row(SIZE, 0.));
		Mat result(SIZE, Row(SIZE, 0.));
		fill_matrix(m1);
		fill_matrix(m2);
		gettimeofday(&setup_end, NULL);
		// setup time
		auto ssec_diff = setup_end.tv_sec * 1. + setup_end.tv_usec * (1e-6) - (setup_start.tv_sec * 1. + setup_start.tv_usec * (1e-6));
		ssec_diff *= 1000.;
		mean_setup_time += ssec_diff;
		setup_times.push_back(ssec_diff);
	}

	mean_setup_time /= iteration * 1.;
	double stand_dev_setup = 0.;
	for (size_t i = 0; i < setup_times.size(); i++)
	{
		stand_dev_setup += (setup_times[i] - mean_setup_time) * (setup_times[i] - mean_setup_time) / setup_times.size() * 1.0;
	}
	stand_dev_setup = sqrt(stand_dev_setup);
	std::cout << "Setup=  " << mean_setup_time << "ms ; Standard Deviation=  " << stand_dev_setup
			  << "ms" << std::endl;

	// matrix multiplication

	// measure performance time for mat mul
	float mean_exec_time = 0;
	// Fills matrices
	Mat m1(SIZE, Row(SIZE, 0.));
	Mat m2(SIZE, Row(SIZE, 0.));
	Mat result(SIZE, Row(SIZE, 0.));
	fill_matrix(m1);
	fill_matrix(m2);
	std::vector<float> exec_times;
	for (int i = 0; i < iteration; i++)
	{
		timeval t1, t2;
		gettimeofday(&t1, NULL);
		mul_matrix(m1, m2, result);
		gettimeofday(&t2, NULL);
		auto sec_diff = t2.tv_sec * 1. + t2.tv_usec * (1e-6) - (t1.tv_sec * 1. + t1.tv_usec * (1e-6));
		sec_diff *= 1000.;
		mean_exec_time += sec_diff;
		exec_times.push_back(sec_diff);
	}

	mean_exec_time /= iteration * 1.;
	double stand_dev = 0.;
	for (size_t i = 0; i < exec_times.size(); i++)
	{
		stand_dev +=
			(exec_times[i] - mean_exec_time) * (exec_times[i] - mean_exec_time) / exec_times.size() * 1.0;
	}
	stand_dev = sqrt(stand_dev);
	std::cout << "Exec=  " << mean_exec_time << "ms ; Standard Deviation=  " << stand_dev
			  << "ms" << std::endl;

	return 0;
}
