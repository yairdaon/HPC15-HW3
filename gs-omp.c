#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

// declare functions
void gsIteration(double * f, double * u ,  int n, double hneg2);
double resNorm(double * u , double * f , int n, double hneg2);
int fullGS(int n);


int main (int argc, char **argv)
{
  

  int n = atoi(argv[1]);
  int gs, t;
  clock_t start, diff;
  printf("\n\n\n");
  printf("Preparing to run. n = %d. \n", n);
  
  start = clock();
  gs  = fullGS(n);
  diff = clock() - start;
  t = (double)diff / CLOCKS_PER_SEC;
  printf("Gauss Seidel. n = %d used. %d iterations. Elapsed time =  %d seconds.\n", n, gs, t);

  return 0;
}




int fullGS(int n) {
 
  printf("\n\n\n");


  // setting arguments
  int gsit = 0;

  double  *f, *gs;
  int i;
  double hneg2 = (n+1) * (n+1) ;
  double initial;
  double factor = 1e6;
  double gsResNorm;

  f  = (double *) malloc(sizeof(double) * n  );
  gs = (double *) malloc(sizeof(double) * n  );



  // set values and stuff
  for( i = 0 ; i < n ; i++) {
    f[i] = 1.0;
    gs[i] = 0.0;
  }

  
  initial = resNorm(gs, f, n, hneg2);
  gsResNorm  = initial;

  // do gauss Seidel iterations, one per round
  while( gsResNorm > initial/factor ) {
      
    gsIteration( f, gs, n, hneg2);
    gsit++;
   
    gsResNorm = resNorm(gs , f, n, hneg2);
    if( gsit % 500000  == 0) 
      printf("Gauss Seidel: ||Au^(%d) - f|| = %f \n", gsit, gsResNorm );
  }  // end while
  
  

  // print results
  printf("Results: \n");
  printf("Gauss Seidel: %d iterations, ||Au - f|| = %f \n", gsit, gsResNorm );

  // free everything
  free(gs);
  free(f);
  return gsit;

}




void gsIteration( double * f, double * u ,  int n, double hneg2){
  /*
    f is an nx1 vector,
    u is an nx1 vector.
  */


  //iteration indices
  int i;
  
  u[0] = (f[0] + hneg2*u[1]) / (2.0*hneg2); 

    
#pragma omp parallel
  {
    
#pragma omp for 
    // update "reds"
    for( i = 1; i < n-1; i+=2){
      u[i] = (f[i] + hneg2*u[i+1] + hneg2*u[i-1]) / (2.0*hneg2);    // old iteration
    }
  
#pragma omp barrier

#pragma omp for 
    // update "blacks"
    for( i = 2; i < n-1; i+=2){
      u[i] = (f[i] + hneg2*u[i+1] + hneg2*u[i-1]) / (2.0*hneg2);    
    }
    
  } // end of parallel

  u[n-1] = (f[n-1] + hneg2*u[n-2]) / (2.0*hneg2);
    
} // end of function scope


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
