/******************************************************************************
* FILE: omp_bug2.c
* DESCRIPTION:
*   Another OpenMP program with a bug. 
* AUTHOR: Blaise Barney 
* LAST REVISED: 04/06/05 
******************************************************************************/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[]) 
{
  int nthreads, i, tid;
  int n = 50;
  float total;

  /*** Spawn parallel region ***/
#pragma omp parallel 
  {
    /* Obtain thread number */
    tid = omp_get_thread_num();
    /* Only master thread does this */
    if (tid == 0) {
      nthreads = omp_get_num_threads();
      printf("Number of threads = %d\n", nthreads);
    }
    printf("Thread %d is starting...\n",tid);

#pragma omp barrier
  
    /* do some work */
    total = 0.0;
#pragma omp parallel for reduction(+:total)
    for (i=0; i < n; ++i) {
      total +=  i*1.0;
    }

#pragma omp barrier
  printf ("done! Total = %e\n",total);  
  printf ("Should be   = %e\n",  (n-1.0)*n/2.0       );
  

  } /*** end parallel ***/

}
// added a reduce command and added barrier before print statement
// for some reason i can't get the thread ids right...
