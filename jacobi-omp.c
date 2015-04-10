#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

// declare function
void jacobiIteration(double *f, double * u, double *v, int n, double hneg2);
double resNorm(double * u , double * f , int n, double hneg2);
int fullJacobi(int n);


int main (int argc, char **argv)
{
 

  int n = 500;//atoi(argv[1]);
  int  jac, t;
  clock_t start, diff;
   
  start = clock();
  jac  = fullJacobi(n);
  diff = clock() - start;
  t = (double)diff / CLOCKS_PER_SEC;
  printf("Jacobi.  n = %d. used %d iterations. Elapsed time =  %d seconds.\n", n, jac, t);

  return 0;
}
  
int fullJacobi(int n){

  printf("\n\n\n");



  // INITIALIZATION STEPS!!!!


  // getting arguments
  int jit  = 0;

  double  *u, *v ,*f;
  int i;
  double hneg2 = (n+1) * (n+1);
  double initial;
  double factor = 1e6;
  double jacResNorm;

  u  = (double *) malloc(sizeof(double) * n  );
  v  = (double *) malloc(sizeof(double) * n  );
  f  = (double *) malloc(sizeof(double) * n  );



  // set values and stuff
  for( i = 0 ; i < n ; i++) {
    f[i] = 1.0;
    u[i] = 0.0;
  }

    
  initial = resNorm(u, f, n, hneg2);
  jacResNorm = initial;


  // FINISHED INITIALIZATION - NOW THE NUMERICS!!!


  // do jacobi iterations, two per round
  while( jacResNorm > initial/factor ) {
   
    jacobiIteration( f, u, v, n, hneg2);
    jit++;
    jacResNorm = resNorm(v, f, n, hneg2);
    //if( jit % 100 == 0 || jit-1 % 100 == 0 )
    //printf("Jacobi: ||Au^(%d) - f|| = %f \n", jit, jacResNorm );

    jacobiIteration(f, v, u, n, hneg2); 
    jit++;
    jacResNorm = resNorm(u, f, n, hneg2);
  }

 


  // print results
  printf("Results: \n");
  printf("Jacobi: %d iterations,  ||Au - f|| = %f \n",  jit, jacResNorm );

  // free everything
  free(f);
  free(v);
  free(u);
  return jit;

}




void jacobiIteration( double * f, double * u, double * v, int n, double hneg2){
  /*
    the current value is at u, the next value is written
    to v. 
    f is an nx1 vector,
    u, v are also nx1
  */

  //iteration indices
  int i; 

  v[0] = (f[0] + hneg2*u[1]) / (2.0*hneg2);


  // the iteration!!!
#pragma omp parallel 
  {

#pragma omp for
  for( i = 1; i < n-1; i++){
    v[i] = (f[i] + u[i-1]*hneg2 + u[i+1]*hneg2) / (2.0*hneg2);
  }

  } // end of parallel part

  v[n-1] = (f[n-1] + hneg2*u[n-2]) / (2.0*hneg2);
 
  
}


double resNorm(double * u, double * f,  int n, double hneg2) {

  int i;
  double norm, tmp;

  tmp  = f[0] - 2.0*hneg2*u[0] + hneg2*u[1];
  norm = tmp*tmp;
#pragma omp parallel
  {

#pragma omp for private(tmp) reduction(+:norm)
  for(i = 1 ; i < n-1 ; i++){
    tmp = f[i] - 2.0*hneg2*u[i] + hneg2*u[i-1] + hneg2*u[i+1];
    norm += tmp*tmp;
  }

  } // end parallel region
  tmp = f[n-1] -2.0*hneg2*u[n-1] + hneg2*u[n-2];
  norm += tmp*tmp;
 
  return sqrt(norm);
}

