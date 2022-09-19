#include <iostream>
#include <cmath>
#include "WLP2D.h"

void cal_Jz_coefficient(double *Jzp){
  double *f = new double[P+1];
  for(int n = 0; n <= M; n++){
    Fai(s*n*dt,f);
    for(int i = 0; i <= P; i++){
      Jzp[i] += Jz(n*dt)*f[i];  
    }
  }
  for(int i = 0; i <= P; i++){
    Jzp[i] = Jzp[i]*s*dt;
  }

  delete [] f;
}
