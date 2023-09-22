#include <iostream>
#include <algorithm>
#include <ctime> 
#include <cstdio>
#include <cmath>

extern "C" {
  #include "blas.h"
}

// Function for convenient printing of matrices:
void print_matrix(const double * const * M, unsigned int m, unsigned int n)
{
    // m is the number of rows, n is the number of columns
    
    // This allows me to avoid typing std:: over and over again.
    // Be careful where you use this though.
    using namespace std;

    // Print out each entry of the matrix:
    cout << "[";
    for(unsigned int i=0; i<m; i++)
    {
      cout << "[";
      for(unsigned int j=0; j<n; j++)
      {
          cout << M[i][j] << ' ';
      }
      cout << "]" << endl;
    }
    cout << "]" << endl;
}

// Useful function for careful deleting of matrices:
void delete_matrix(double **a, unsigned int m)
{
    //m is the number of rows of a
    for(unsigned int i=0; i<m; i++)
    {
        delete [](a[i]);
    }
    delete []a;
}

// This function returns a random matrix of size m x n
double** generate_random_matrix(unsigned int m, unsigned int n)
{
    // m is the number of rows, n is the number of columns
    
    // This allows me to avoid typing std:: over and over again.
    // Be careful where you use this though.
    using namespace std;

    //Declare the matrix:
    double** M = new double*[m];

    // rand produces a random INTEGER between 0 and RAND_MAX.
    //   We have to divide by RAND_MAX to get a random double between 0 and 1
    // Note that curly braces are not necessary for 1-line for-loops but
    //   I've left them in to avoid confusion.
    for(unsigned int i=0; i<m; i++)
    {
        M[i] = new double[n];
        for(unsigned int j=0; j<n; j++)
        {
             M[i][j] = ((double) rand()/(RAND_MAX));
        }
    }

    return M;
}

// This function returns a zero matrix of size m x n
double** generate_zero_matrix(unsigned int m, unsigned int n)
{
    // m is the number of rows, n is the number of columns
    
    // This allows me to avoid typing std:: over and over again.
    // Be careful where you use this though.
    using namespace std;

    //Declare the matrix:
    double** M = new double*[m];

    // rand produces a random INTEGER between 0 and RAND_MAX.
    //   We have to divide by RAND_MAX to get a random double between 0 and 1
    // Note that curly braces are not necessary for 1-line for-loops but
    //   I've left them in to avoid confusion.
    for(unsigned int i=0; i<m; i++)
    {
        M[i] = new double[n];
        for(unsigned int j=0; j<n; j++)
        {
             M[i][j] = 0.0;
        }
    }

    return M;
}

// Compare two matrices:
bool compare_matrices(const double * const * M1, const double * const * M2,
                      unsigned int m, unsigned int n, double tol = 1.0e-11)
{
    // m is the number of rows, n is the number of columns
    for(unsigned int i=0; i<m; i++)
        for(unsigned int j=0; j<n; j++)
            if( std::abs( M1[i][j]-M2[i][j] ) > tol )
              return false;
    return true;
}

// dgemm only takes in 1D arrays so we need this to flatten it
double* flatten_matrix(const double * const * M, unsigned int m, unsigned int n)
{
    double *M_flat = new double[m*n];
    //Note that fortran is column-major
    for(unsigned int i=0; i<n; i++)
        for(unsigned int j=0; j<m; j++)
             M_flat[i*m+j] = M[j][i];

    return M_flat;
}

void call_dgemm(int m, int n, int k,
              double alpha, const double * const * A, 
              const double * const * B, 
              double beta, double **C )
{
    double *A_flat = flatten_matrix(A, m, k);
    double *B_flat = flatten_matrix(B, k, m);
    double *C_flat = flatten_matrix(C, m, n);

    char matrixType = 'n';
/*
    dgemm_(&matrixType, &matrixType, &m, &n, &k, 
            &alpha, A_flat, &m, B_flat, &k,
            &beta, C_flat, &m);
*/
    //Unflatten solution:
    for(unsigned int i=0; i<n; i++)
        for(unsigned int j=0; j<m; j++)
             C[j][i] = C_flat[i*m+j];

    delete []A_flat;
    delete []B_flat;
    delete []C_flat;
}


// Fill in this function with a matrix-matrix multiply routine:
void mydgemm( char transa, char transb,
              unsigned int m, unsigned int n, unsigned int k, 
              double alpha, const double * const * a,
              const double * const * b,
              double beta, double **c)
{
    //This function should perform the following operation:
    //C = alpha*A*B + beta*C

    /*
     Simplifying assumptions
        * Assume that matrices a, b, and c are all square (m == n == k)
          This also means lda == ldb == ldc == n
        * Ignore transa and transb for this exercise. Just pass in 'n'
          for these arguments.
        * You do not need to perform multiplication by alpha or beta.
          You may assume these are 1 and 0, respectively.

     Inputs:
      transa = assume 'n'
      transb = assume 'n'
           m = # rows of a and c
           n = # columns in b and c
           k = # columns of a and rows of b
       alpha = assume 1.0
           a = the matrix A
         lda = leading dimension of A
           b = the matrix B
         ldb = leading dimension of B
           c = the matrix C
         ldc = leading dimension of C

     Output:
      An updated matrix c
    */

    //Write your naive matrix multiply here
  
}

int main(int argc, char** argv)
{
    //Size of matrix
    int n;
    n=0;
    if( argc > 1)
    {
        n = atoi(argv[1]);
    }
    if( n <= 0 )
    {
        std::cout << n << " is an invalid size!" << std::endl;
        return 1;
    }
    std::cout << " Performing matrix-matrix multipliciation for SIZE=" << n << std::endl;

    double **a = generate_random_matrix(n,n);
    double **b = generate_random_matrix(n,n);
    double **c = generate_zero_matrix(n,n);
    double **d = generate_zero_matrix(n,n);

    // This is only a valid command because of our overloading of << in
    //   the templated function above.
    //std::cout << "a = ";
    //print_matrix(a,n,n);
    //std::cout << "b = ";
    //print_matrix(b,n,n);
    //std::cout << "c = ";
    //print_matrix(c,n,n);
    //std::cout << "d = ";
    //print_matrix(d,n,n);

    clock_t start;
    double duration;

    //Record the current time:
    start = clock(); 
    // Perform c = a*b+c:
    mydgemm( 'n', 'n', n, n, n, 1.0, a, b, 1.0, c );
    // Compute how long it took:
    duration = ( clock()-start ) / (double)(CLOCKS_PER_SEC);

    //std::cout << std::endl << "a*b = ";
    //print_matrix(c,n,n);
    std::cout << std::endl << " mydgemm took " << duration << 
                              " seconds!" << std::endl;

//  Uncomment this to call blas:
//    //Call blas:
//    start = clock(); 
//    call_dgemm(n,n,n,1.0,a,b,1.0,d);
//    duration = ( clock()-start ) / (double)(CLOCKS_PER_SEC);
//    //std::cout << std::endl << "a*b = ";
//    //print_matrix(d,n,n);
//    std::cout << std::endl << "BLAS dgemm took " << duration << 
//                              " seconds!" << std::endl;
//
//    //Check to see if answer is correct:
//    if(!compare_matrices(c,d,n,n))
//      std::cout << "Uh oh! Is your matrix multiply correct?" << std::endl;


    delete_matrix(a,n);
    delete_matrix(b,n);
    delete_matrix(c,n);
    delete_matrix(d,n);
    return 0;
}
