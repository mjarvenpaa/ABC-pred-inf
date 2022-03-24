#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

// Simulates one step of the 3-parameter version of LV
void LV3_simul_step(int *xy, double *t0p, double *dtp, double *th)
{
  double t=*t0p, dt=*dtp, tend=t+dt;
  double hsum,h1,h2,h3,u;
  GetRNGstate();
  while (1) {
    h1=th[0]*xy[0]; h2=th[1]*xy[0]*xy[1]; h3=th[2]*xy[1];
    hsum=h1+h2+h3;
    if ((hsum<1e-6)||(xy[0]>1e5)) {
      PutRNGstate();
      return;
    } else {
      t+=rexp(1.0/hsum);
    }
    if (t>=tend) {
      PutRNGstate();
      return;
    }
    u=unif_rand();
    if (u<h1/hsum) {
      xy[0]+=1;
    } else if (u<(h1+h2)/hsum) {
      xy[0]-=1;
      xy[1]+=1;
    } else {
      xy[1]-=1;
    }
  }
}

// Simulates one step of the 4-parameter version of LV
void LV4_simul_step(int *xy, double *t0p, double *dtp, double *th)
{
  double t=*t0p, dt=*dtp, tend=t+dt;
  double hsum,h1,h2,h3,h4,u;
  GetRNGstate();
  while (1) {
    h1=th[0]*xy[0]; h2=th[1]*xy[0]*xy[1]; h3=th[2]*xy[0]*xy[1]; h4=th[3]*xy[1];
    hsum=h1+h2+h3+h4;
    if ((hsum<1e-6)||(xy[0]>1e5)) {
      PutRNGstate();
      return;
    } else {
      t+=rexp(1.0/hsum);
    }
    if (t>=tend) {
      PutRNGstate();
      return;
    }
    u=unif_rand();
    if (u<h1/hsum) {
      xy[0]+=1;
    } else if (u<(h1+h2)/hsum) {
      xy[0]-=1;
    } else if (u<(h1+h2+h3)/hsum) {
      xy[1]+=1;
    } else {
      xy[1]-=1;
    }
  }
}

/* Full simulation of LV model, implements the following R loop:
 xy <- matrix(0,2,nn)
 xy[,1] <- xy0
 for (i in 2:nn) {
  if (np==3) {
    xy[,i] <- LV3.simul.step(xy[,i-1], tobs[i-1], tobs[i] - tobs[i-1], th)
  } else if (np==4) {
    xy[,i] <- LV4.simul.step(xy[,i-1], tobs[i-1], tobs[i] - tobs[i-1], th)
  } else {
    stop('Incorrect parameter for LV model.')
  }
 }
 The implementation is a bit clumsy but this allows e.g. LV3_simul_step to be
 also directly called from R. Using .Call(...) or Rcpp would allow less clumsy 
 implementation but this now suffices here.
*/
void LV_simul_full(int *x, int *y, int *xy0, double *tps, int *nn, double *th, int *np)
{
  int i;
  double dti,t0pi;
  int xy[2];
  xy[0] = xy0[0]; xy[1] = xy0[1]; 
  x[0] = xy0[0]; y[0] = xy0[1]; 
  for (i=1; i<*nn; ++i) {
    t0pi=tps[i-1];
    dti=tps[i]-tps[i-1];
    if (*np==3) {
      LV3_simul_step(xy, &t0pi, &dti, th);
    } else if (*np==4) {
      LV4_simul_step(xy, &t0pi, &dti, th);
    } else {
      error("Incorrect parameter for LV model.");
    }
    x[i] = xy[0]; y[i] = xy[1]; 
  }
}


