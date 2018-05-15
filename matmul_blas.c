/*
 * Université Pierre et Marie Curie
 *
 * Programme de multiplication de matrices carrees.
 */ 

#include <stdlib.h>
#include <stdio.h>

#include <sys/time.h>

#include <common.h>
/* common.h inclut common_x86_64.h 
   Pour le vÃ©rifier : 
   #ifdef ARCH_X86_64
   #error "toto"
   #endif
*/
#include <cblas.h>

double my_gettimeofday(){
  struct timeval tmp_time;
  gettimeofday(&tmp_time, NULL);
  return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6L);
}


#define REAL_T float 
#define NB_TIMES 10

/*** Matmul: ***/
/* C += A x B 
 * square matrices of order 'n'
 */
void matmul(int n, REAL_T *A, REAL_T *B, REAL_T *C){
  /* CBLAS interface : C = alpha A*B + beta*C 
     void cblas_sgemm(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, 
     blasint M, blasint N, blasint K,
     float alpha, float *A, blasint lda, 
     float *B, blasint ldb, 
     float beta, float *C, blasint ldc);
  */
  cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
	      n, n, n, 
	      1.0 /* alpha */, A, n, 
	      B, n, 
	      1.0 /* beta */, C, n);
}


int main(int argc, char **argv)
{
  int i,j;
  double debut=0.0, fin=0.0;
  REAL_T *A, *B, *C;
  int n=2; /* default value */
  int nb=0;
  
  /* Read 'n' on command line: */
  if (argc == 2){
    n = atoi(argv[1]);
  }

  /* Allocate the matrices: */
  if ((A = (REAL_T *) malloc(n*n*sizeof(REAL_T))) == NULL){
    fprintf(stderr, "Error while allocating A.\n");
  }
  if ((B = (REAL_T *) malloc(n*n*sizeof(REAL_T))) == NULL){
    fprintf(stderr, "Error while allocating B.\n");
  }
  if ((C = (REAL_T *) malloc(n*n*sizeof(REAL_T))) == NULL){
    fprintf(stderr, "Error while allocating C.\n");
  }

  /* Initialize the matrices */
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++){
      *(A+i*n+j) = 1 / ((REAL_T) (i+j+1));
      *(B+i*n+j) = 1.0;
      *(C+i*n+j) = 1.0;
    }

  /* Start timing */
  debut = my_gettimeofday();
  for (nb=0; nb<NB_TIMES; nb++){
    /* Do matrix-product C=A*B+C */
    matmul(n, A, B, C);
    /* End timing */
  }
  fin = my_gettimeofday();

  fprintf( stdout, "For n=%d: total computation time (with gettimeofday()) : %g s\n",
	   n, (fin - debut)/NB_TIMES);
  fprintf( stdout, "For n=%d: performance = %g Gflop/s \n",
	   n, (((double) 2)*n*n*n / ((fin - debut)/NB_TIMES) )/ ((double) 1e9) ); /* 2n^3 flops */
      
  /* Print 2x2 top-left square of C : */
  for(i=0; i<2 ; i++){
    for(j=0; j<2 ; j++)
      printf("%+e  ", C[i*n+j]);
    printf("\n");
  }
  printf("\n");
  /* Print 2x2 bottom-right square of C : */
  for(i=n-2; i<n ; i++){
    for(j=n-2; j<n ; j++)
      printf("%+e  ", C[i*n+j]);
    printf("\n");
  }

  /* Free the matrices: */
  free(A); 
  free(B); 
  free(C);

  return 0;
}
