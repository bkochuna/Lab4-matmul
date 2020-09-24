PROGRAM matmul_example
  IMPLICIT NONE

  !Variable for processing input
  CHARACTER(LEN=4) :: nstr
  INTEGER :: n,ioerr

  !Variables for timing
  INTEGER :: start,finish,count_rate
  REAL(8) :: elapsed_time

  !Variables for matrix multiplication
  REAL(8),ALLOCATABLE :: a(:,:),b(:,:),c(:,:),d(:,:)


  !Get the size of matrix from command line
  nstr=''
  CALL GET_COMMAND_ARGUMENT(1,VALUE=nstr)

  !Default size
  n=0
  
  !Convert to an integer
  READ(nstr,*,IOSTAT=ioerr) n

  !Check your input
  IF(ioerr /= 0) THEN
    WRITE(*,*) 'Command argument "'//TRIM(nstr)// &
      '" does not appear to be an integer!'
    STOP 1
  ENDIF

  IF(n > 0) THEN
    WRITE(*,*) "Performing matrix-matrix multipliciation for SIZE=",n
    WRITE(*,*)

    !Create storage for matrices
    ALLOCATE(a(n,n))
    ALLOCATE(b(n,n))
    ALLOCATE(c(n,n))
    ALLOCATE(d(n,n))
  
    !Fill matrices with random numbers
    CALL RANDOM_NUMBER(a)
    CALL RANDOM_NUMBER(b)
  
    CALL SYSTEM_CLOCK(start,count_rate) !get start time

    !my general matrix-matrix multiply
    CALL mydgemm('n','n',n,n,n,1.0d0,a,n,b,n,0.0d0,c,n)

    CALL SYSTEM_CLOCK(finish) !get finish time

    !Convert time to seconds and print
    elapsed_time=REAL(finish-start,8)/REAL(count_rate,8)
    WRITE(*,'(a,f6.3,a)') " mydgemm took",elapsed_time," seconds"

!    CALL SYSTEM_CLOCK(start,count_rate) !get start time
!
!    !my general matrix-matrix multiply
!    CALL dgemm('n','n',n,n,n,1.0d0,a,n,b,n,0.0d0,d,n)
!
!    CALL SYSTEM_CLOCK(finish) !get finish time
!
!    !Convert time to seconds and print
!    elapsed_time=REAL(finish-start,8)/REAL(count_rate,8)
!    WRITE(*,'(a,f6.3,a)') "BLAS dgemm took",elapsed_time," seconds!"

     !Compare the result of your matrix multiply with BLAS
!    IF(ANY(ABS(c-d) > 1.0d-13)) THEN
!      WRITE(*,*) "Uh oh! Is your matrix multiply correct?"
!      WRITE(*,*) MAXVAL(ABS(c-d))
!    ENDIF
  ELSE
    WRITE(*,*) "Matrix size must be greater than 0!"
    STOP 1
  ENDIF

CONTAINS

  SUBROUTINE mydgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
    CHARACTER(LEN=*),INTENT(IN) :: transa,transb
    INTEGER,INTENT(IN) :: m,n,k,lda,ldb,ldc
    REAL(8),INTENT(IN) :: alpha,beta
    REAL(8),INTENT(IN) :: a(lda,*),b(ldb,*)
    REAL(8),INTENT(INOUT) :: c(ldc,*)

    !This function should perform the following operation:
    !C = alpha*A*B + beta*C

    !Simplifying assumptions
    !   * Assume that matrices a, b, and c are all square (m == n == k)
    !     This also means lda == ldb == ldc == n
    !   * Ignore transa and transb for this exercise. Just pass in 'n'
    !     for these arguments.
    !   * You do not need to perform multiplication by alpha or beta.
    !     You may assume these are 1 and 0, respectively.

    !Inputs:
    ! transa = assume 'n'
    ! transb = assume 'n'
    !      m = # rows of a and c
    !      n = # columns in b and c
    !      k = # columns of a and rows of b
    !  alpha = assume 1.0
    !      a = the matrix A
    !    lda = leading dimension of A
    !      b = the matrix B
    !    ldb = leading dimension of B
    !      c = the matrix C
    !    ldc = leading dimension of C

    !Output:
    !  An updated matrix c

    !Write your naive matrix multiply here

  ENDSUBROUTINE mydgemm

ENDPROGRAM matmul_example
