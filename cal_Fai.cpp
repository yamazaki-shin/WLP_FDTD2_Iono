#include <iostream>
#include <cmath>
#include "WLP2D.h"

void Fai(double t, double *F){
double *L = new double [P+1];

  L[0] = 1.0;
  L[1] = 1.0 - t;
  F[0] = exp(-t/2.0)*L[0];
  F[1] = exp(-t/2.0)*L[1];
  for(int n = 2; n < P+1; n++){
    L[n] = ((2.0*n - 1.0 - t) * L[n-1] - (n-1)*L[n-2])/n;
    F[n] = exp(-t/2.0)*L[n];
  }
  delete [] L;
}
