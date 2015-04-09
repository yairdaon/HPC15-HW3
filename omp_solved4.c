/******************************************************************************
* FILE: omp_bug4.c
* DESCRIPTION:
*   This very simple program causes a segmentation fault.
* AUTHOR: Blaise Barney  01/09/04
* LAST REVISED: 04/06/05
******************************************************************************/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#define N 1048


int main (int argc, char *argv[]) 
{
  printf("\n\nOutput for mpi_solved4.\n\n");

  int nthreads, tid, i, j;
  double a[N*N];

  /* Fork a team of threads with explicit variable scoping */
#pragma omp parallel shared(nthreads) private(i,j,tid)
  {


    /* Obtain/print thread info */
    tid = omp_get_thread_num();
    if (tid == 0) 
      {
	;
	nthreads = omp_get_num_threads();
	printf("Number of threads = %d\n", nthreads);
      }
    printf("Thread %d starting...\n", tid);

    /* Each thread works on its own private copy of the array */
    for (i=0; i<N; i++)
      for (j=0; j<N; j++)
	a[i*N + j] = tid + i + j;
      
    /* For confirmation */
    printf("Thread %d done. Last element= %f\n",tid, a[ N*N-1 ]);

  }  /* All threads join master thread and disband */
  
  printf("\n\n\n");
  return 0;
}

// a is a pointer so just because it is private we didn't really do anything
// Solution is to remove it from the private clause. But this makes no sense!!!
// OK, so it does, kind of. If we declare a to be private we make the *pointer*
// private. Without it, we make the entire array private...?
