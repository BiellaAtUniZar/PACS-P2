#include <assert.h>
#include <chrono>
#include <iostream>
#include <random>
#include <vector>
#include <sys/time.h>
#include <array>
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
 * Fills the square matrix of dimension SIZExSIZE of random values between -10
 * and 10 included
 * @param size matrix size
 */

float* fill_matrix_dynamic(int size)
{
	random_device rd; // Setup of number generator
	mt19937 gen(rd());
	uniform_real_distribution<> double_dist(-10.0, 10.0);
	//float* array = new float[size*size];  
	int k = 0;
	float* array = new float[size*size];
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			// matrix[i][j] = double_dist(gen);  // Filling the array with random
			// doubles
			k++;
		    array[i * size + j] = k;
		}
	}
	return array;
}

/**
 * Performs multiplication and returns address of resulting matrix
 * @param matrix1 
 * @param matrix2 
 * @param size  matrices dimension
 */
float* mul_matrix_dynamic(const float *matrix1, const float *matrix2, int size)
{
    // Allocate result matrix and zero-initialize it
    float* result = new float[size * size]();

    // Perform matrix multiplication
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            float sum = 0.0f;
            for (int k = 0; k < size; k++) {
                // Multiply and accumulate for result[i, j]
                sum += matrix1[i * size + k] * matrix2[k * size + j];
            }
            result[i * size + j] = sum; // Store result
        }
    }

    return result; // Return dynamically allocated result matrix
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

/**
 * Prints  matrix
 * @param matrix1 
 * @param size 
 */

void print_matrix_dynamic(const float *matrix,int size)
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j <size; j++)
		{
			std::cout << matrix[i*size+j] << " ";
		}
		std::cout << std::endl;
	}
}

/**
 * Prints  matrix
 * @param matrix 
 */

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
	bool AUTO_MODE = argc != 2;
	if (AUTO_MODE)
	{
		SIZE = 10;
	}
	else
	{
		SIZE = stoi(argv[1]); // Parses parameter, throws exception if not int
	}
    
	// matrix multiplication
	timeval setup_start, setup_end;
	gettimeofday(&setup_start, NULL);
	float * m1= fill_matrix_dynamic(SIZE);
	float * m2= fill_matrix_dynamic(SIZE);
	gettimeofday(&setup_end, NULL);

	timeval t1, t2;
	gettimeofday(&t1, NULL);
    float *result= mul_matrix_dynamic(m1, m2,SIZE);
	gettimeofday(&t2, NULL);
	
	delete[] m1;
	delete[] m2;
	delete[] result;
	return 0;
}
