/******************************************************************************
* FILE: omp_bug6.c
* DESCRIPTION:
*   This program compiles and runs fine, but produces the wrong result.
*   Compare to omp_orphan.c.
* AUTHOR: Blaise Barney  6/05
* LAST REVISED: 06/30/05
******************************************************************************/

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#define VECLEN 100


float a[VECLEN], b[VECLEN];
float sum = 0.0;

float dotprod ()
{
  int i,tid;

  tid = omp_get_thread_num();
#pragma omp for reduction(+:sum)
  for (i=0; i < VECLEN; i++)
    {
      sum = sum + (a[i]*b[i]);
      //printf("tid= %d i=%d\n",tid,i);
    }
}


int main (int argc, char *argv[]) {
  printf("\n\nOutput for mpi_solved6.\n\n");

  int i;
 

  for (i=0; i < VECLEN; i++)
    a[i] = b[i] = 1.0 * i;
  sum = 0.0;

#pragma omp parallel 
  dotprod();

  printf("Sum = %f. Should be %d\n",sum, VECLEN*(VECLEN+1)*(2*VECLEN+1)/6    );

  printf("\n\n\n");
  return 0;

}

// Per the advice in omp_orphan.c, I mad the sum variable global (declare outside any function).

