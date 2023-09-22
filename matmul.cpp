#include <algorithm>
#include <chrono>
#include <iostream>
#include <random>

extern "C" {
  #include "blas.h"
}

// This function returns a random matrix of size m x n
// m is the number of rows, n is the number of columns
double * generate_random_matrix(int m, int n)
{

  // Declare the matrix:
  // We want to allocate a single block of memory of size m * n
  // instead of allocating m blocks of size n. This is called
  // a "flattened" matrix, and it ensures that the entire matrix
  // is contiguous in memory. This can help with performance.
  double * matrix = new double[m * n];

  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_real_distribution<double> dist(0.0, 1.0);

  for (int i = 0; i < m * n; ++i) {
    matrix[i] = dist(rng);
  }
  // or for a more modern approach:
  // std::generate(matrix, matrix + m * n, [&]() { return dist(rng); });

  return matrix;
}

// This function returns a zero matrix of size m x n
// m is the number of rows, n is the number of columns
double * generate_zero_matrix(int m, int n)
{
  double * matrix = new double[m * n];

  for (int i = 0; i < m * n; ++i) {
    matrix[i] = 0.0;
  }
  // or for a more modern approach:
  // std::fill(matrix, matrix + m * n, 0.0);

  return matrix;
}

// Fill in this function with a matrix-matrix multiply routine:
void mydgemm(char transa, char transb,
              int m, int n, int k,
              double alpha,
              double const * const a,
              double const * const b,
              double beta,
              double * c)
{
  // This function should perform the following operation:
  // C = alpha*A*B + beta*C

  // Simplifying assumptions
  //  * Assume that matrices a, b, and c are all square (m == n == k)
  //    This also means lda == ldb == ldc == n
  //  * Ignore transa and transb for this exercise. Just pass in 'n'
  //    for these arguments.
  //  * You do not need to perform multiplication by alpha or beta.
  //    You may assume these are 1 and 0, respectively.

  // Inputs:
  //  transa = assume 'n'
  //  transb = assume 'n'
  //       m = # rows of a and c
  //       n = # columns in b and c
  //       k = # columns of a and rows of b
  //   alpha = assume 1.0
  //       a = the matrix A
  //     lda = leading dimension of A
  //       b = the matrix B
  //     ldb = leading dimension of B
  //       c = the matrix C
  //     ldc = leading dimension of C

  // Output:
  //  An updated matrix c

  // Write your naive matrix multiply here

}

void transpose_matrix(double * const matrix, int m, int n)
{
  // This function transposes a matrix in place.
  // The matrix is m x n in size.
  for (int i = 0; i < m; ++i) {
    for (int j = i + 1; j < n; ++j) {
      std::swap(matrix[i * n + j], matrix[j * n + i]);
    }
  }
}

void call_dgemm(char transa, char transb,
              int m, int n, int k,
              double alpha,
              double * const a,
              double * const b,
              double beta,
              double * const c)
{
  dgemm_(&transa, &transb, &m, &n, &k,
          &alpha, a, &m, b, &k, &beta, c, &m);
}

// Compare two matrices:
bool compare_matrices(double const * const m1, double const * const m2,
                      int m, int n, double tol = 1.0e-11)
{
  for (int i = 0; i < m * n; ++i) {
    if (std::abs(m1[i] - m2[i]) > tol) {
      return false;
    }
  }
  return true;
}

int main(int argc, char** argv)
{
  //  Size of tje matrix
  int n = 0;
  if (argc > 1) {
    n = atoi(argv[1]);
  }
  if( n <= 0 ) {
    std::cerr << n << " is an invalid size!" << std::endl;
    return 1;
  }
  std::cout << " Performing matrix-matrix multipliciation for SIZE=" << n << std::endl;

  double * a = generate_random_matrix(n, n);
  double * b = generate_random_matrix(n, n);
  double * c = generate_zero_matrix(n, n);
  double * d = generate_zero_matrix(n, n);

  //  Record the current time:
  auto start = std::chrono::high_resolution_clock::now();

  // Perform C = A * B
  mydgemm('n', 'n', n, n, n, 1.0, a, b, 1.0, c );

  // Compute how long it took:
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration = end - start;
  std::cout <<"\n   mydgemm took " << duration.count() << " seconds!" << std::endl;

  // Due to row-major ordering of C/C++ and column-major ordering of Fortran,
  // we need to transpose the inputs when calling blas to get the same result.
  // In our dgemm, we compute C = A * B.
  // But, if we note that C^T = (A * B)^T = B^T * A^T, then we can see that
  // we just need to swap the order of the inputs to get the same result
  start = std::chrono::high_resolution_clock::now();
  call_dgemm('n', 'n', n, n, n, 1.0, b, a, 1.0, d);
  // Or we can transpose the inputs and call dgemm like this:
  // call_dgemm('t', 't', n, n, n, 1.0, a, b, 1.0, d);
  // But the result will be transposed, so we need to transpose it again
  end = std::chrono::high_resolution_clock::now();
  duration = end - start;
  std::cout <<"BLAS dgemm took " << duration.count() << " seconds!" << std::endl;

  //Check to see if answer is correct:
  if(!compare_matrices(c, d, n, n)) {
    std::cout << " Uh oh! Is your matrix multiply correct?" << std::endl;
  }

  delete[] a;
  delete[] b;
  delete[] c;
  delete[] d;
  return 0;
}
