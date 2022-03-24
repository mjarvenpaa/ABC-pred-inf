#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

/* Implements the following R loop of the Simple Markov model simulation:
 x[1] <- cp + es[1]
 if (n>1) {
  for (i in 2:n) {
    x[i] = cp + phi*x[i-1] + es[i]
  }
 }
*/
void simplemarkov_simul_loop(double *x, double *es, double *cp, double *phi, int *n)
{
  int i;
  x[0] = *cp + es[0];
  for (i=1; i<*n; ++i) {
    x[i] = *cp + *phi*x[i-1] + es[i];
  }
}
