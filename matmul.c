#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "blas.h"

// This function returns a random matrix of size m x n
// m is the number of rows, n is the number of columns
double ** generate_random_matrix(int m, int n)
{

  // Declare the matrix:
  // Note that this is not a contiguous allocation! This will 
  // have less efficient memory access patterns than a contiguous
  // allocation.
  double ** matrix = (double **)malloc(m * sizeof(double *)); 
  
  for (int i = 0; i < m; ++i) {
    matrix[i] = (double *)malloc(n * sizeof(double));
    for (int j = 0; j < n; ++j) {
      matrix[i][j] = (double)rand() / (double)RAND_MAX; 
    }
  }

  return matrix;
}

// This function returns a zero matrix of size m x n
// m is the number of rows, n is the number of columns
double ** generate_zero_matrix(int m, int n)
{
  double ** matrix = (double **)malloc(m * sizeof(double *)); 
  for (int i = 0; i < m; ++i) {
    // We can use calloc here to zero-initialize the matrix
    matrix[i] = (double *)calloc(n, sizeof(double));
  }

  return matrix;
}

// Fill in this function with a matrix-matrix multiply routine:
void mydgemm(char transa, char transb,
              int m, int n, int k,
              double alpha,
              double ** a,
              double ** b,
              double beta,
              double ** c)
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

double * flatten_matrix(double ** matrix, int m, int n)
{
  // This function takes a 2D matrix and returns a 1D array
  // with the same data. This is useful for passing the matrix
  // to BLAS functions.
  double * flat = (double *)malloc(m * n * sizeof(double));
  for (int i = 0; i < m; ++i) {
    memcpy(flat + i * n, matrix[i], n * sizeof(double));
  }
  return flat;
}

void unflatten_matrix(double * flat, double ** matrix, int m, int n)
{
  for (int i = 0; i < m; ++i) {
    memcpy(matrix[i], flat + i * n, n * sizeof(double));
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
int compare_matrices(double ** m1, double ** m2,
                      int m, int n, double tol)
{
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; j++) {
      if (abs(m1[i][j] - m2[i][j]) > tol) {
        return 0; 
      }
    }
  }
  return 1;
}

int main(int argc, char** argv)
{
  //  Size of tje matrix
  int n = 0;
  if (argc > 1) {
    n = atoi(argv[1]);
  }
  if (n <= 0) {
    fprintf(stderr, "%s is an invalid size!\n", argv[1]);
    return 1;
  }
  printf(" Performing matrix-matrix multiplication for SIZE=%d\n", n);

  double ** a = generate_random_matrix(n, n);
  double ** b = generate_random_matrix(n, n);
  double ** c = generate_zero_matrix(n, n);
  double ** d = generate_zero_matrix(n, n);

  //  Record the current time:
  clock_t start = clock();

  // Perform C = A * B
  mydgemm('n', 'n', n, n, n, 1.0, a, b, 1.0, c );

  // Compute how long it took:
  clock_t end = clock();
  double elapsed_time = (end - start) / (double)CLOCKS_PER_SEC;
  printf("\n   mydgemm took %f seconds!\n", elapsed_time);

  // Due to row-major ordering of C/C++ and column-major ordering of Fortran,
  // we need to transpose the inputs when calling blas to get the same result.
  // In our dgemm, we compute C = A * B.
  // But, if we note that C^T = (A * B)^T = B^T * A^T, then we can see that
  // we just need to swap the order of the inputs to get the same result
  //
  // BLAS expects flat arrays, so:
  double * a_flat = flatten_matrix(a, n, n);
  double * b_flat = flatten_matrix(b, n, n);
  double * d_flat = flatten_matrix(d, n, n);
  start = clock();
  call_dgemm('n', 'n', n, n, n, 1.0, b_flat, a_flat, 1.0, d_flat);
  // Or we can transpose the inputs and call dgemm like this:
  // call_dgemm('t', 't', n, n, n, 1.0, a, b, 1.0, d);
  // But the result will be transposed, so we need to transpose it again
  end = clock(); 
  unflatten_matrix(d_flat, d, n, n);
  elapsed_time = (end - start) / (double)CLOCKS_PER_SEC;
  printf("BLAS dgemm took %f seconds!\n", elapsed_time);

  //Check to see if answer is correct:
  if(!compare_matrices(c, d, n, n, 1.0e-11)) {
    printf(" Uh oh! Is your matrix multiply correct?\n");
  }

  free(a_flat);
  free(b_flat);
  free(d_flat);
  free(a);
  free(b);
  free(c);
  free(d);
  return 0;
}
