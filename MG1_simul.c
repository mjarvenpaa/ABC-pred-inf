#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

/* Implements the following loop of the M/G/1 model simulation:
 if (n>1) {
    for (i in 2:n) {
      y[i] <- u[i] + max(0, v[i] - x[i-1])
      x[i] <- x[i-1] + y[i]
    }
 }
 The beginning of M/G/1 simulation, i.e. pre-simulation of the random variables 
 w etc. and computing the first element in x and y, is not implemented here for 
 simplicity and because it would likely not produce much improvement as it is
 already fast/vectorized in R.
 */
void MG1_simul_loop(double *x, double *y, double *v, double *u, int *n)
{
  int i;
  for (i=1; i<*n; ++i) {
    y[i] = u[i] + fmax(0.0, v[i] - x[i-1]);
    x[i] = x[i-1] + y[i];
  }
}
