void dgemm_(
    char * transa, char * transb,
    int * m, int * n, int * k,
    double * alpha,
    double * A, int * lda,
    double * B, int * ldb,
    double * beta,
    double * C, int * ldc);
