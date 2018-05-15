/*
 * Université Pierre et Marie Curie
 *
 * Programme de multiplication de matrices carrees.
 */ 

#include <stdlib.h>
#include <stdio.h>

#include <sys/time.h>

double my_gettimeofday(){
  struct timeval tmp_time;
  gettimeofday(&tmp_time, NULL);
  return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6L);
}


#define REAL_T float 
#define NB_TIMES 10

/* Data L1 cache size = 32 Ko 
 * We only consider here B (not the sizes of the 3 matrices), 
 * which implies  8196 floats => threshold = 8192   
 * (sqrt(8192) = 90.5 -> donc ok pour mes valeurs de N testees) 
 *
 * L2 cache size : 2 Mo partagÃ©s par 2 coeurs 
 * 512 * 1024 floats => threshold = 524288
 * (sqrt(512*1024) = 724.1 -> donc ok pour mes valeurs de N testees) */
int threshold = 8192; // -> perf ok jusqu'a 1000
//int threshold = 524288; // -> perf pas bonnes...

/*** 'mm' recursive function : ***/ 
void mm (int crow, int ccol, /* Corner of C block */
	 int arow, int acol, /* Corner of A block */
	 int brow, int bcol, /* Corner of B block */
	 int n,              /* Blocks are n x n */
	 int stride,         /* Stride for whole matrix (original value of 'n') */
	 REAL_T *A, REAL_T *B, REAL_T *C){
  int nhalf[3]; /* Quadrant sizes */
  int i,j,k; 
  REAL_T *aptr, *bptr, *cptr;
  
  if (n * n > threshold){
    /* B does not fit in cache -> multiply blocks of A, B */
    
    nhalf[0] = 0; nhalf[1] = n/2 ; nhalf[2] = n - n/2; 
    
    for (i=0; i<2; i++){
      for (j=0; j<2; j++){
	for (k=0; k<2; k++){
	  mm(crow+nhalf[i], ccol+nhalf[j], 
	     arow+nhalf[i], acol+nhalf[k],
	     brow+nhalf[k], bcol+nhalf[j],
	     nhalf[j+1], stride,
	     A, B, C);
	}
      }
    }

  }
  else {
    /* B fits in cache -> do standard multiply */
    
    for (i=0; i<n; i++){
      for (j=0; j<n; j++){
	cptr = &C[(crow+i)*stride + ccol+j];
	aptr = &A[(arow+i)*stride + acol];
	bptr = &B[(brow  )*stride + bcol+j];
	for (k=0; k<n; k++){
	  *cptr += *(aptr++) * *bptr; 
	  bptr += stride;
	}
      }
    }

  }
}



/*** Matmul: ***/
/* C += A x B 
 * square matrices of order 'n'
 */
void matmul(int n, REAL_T *A, REAL_T *B, REAL_T *C){

  mm(0, 0, 
     0, 0, 
     0, 0, 
     n, n /* stride */,
     A, B, C);
}


int main(int argc, char **argv)
{
  int i,j;
  double debut=0.0, fin=0.0;
  REAL_T *A, *B, *C;
  int n=2; /* default value */
  int nb=0;
  
  /* Read 'n' and 'threshold' on command line: */
  if (argc >= 2){
    n = atoi(argv[1]);
  }
  if (argc >= 3){
    threshold = atoi(argv[2]);
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
