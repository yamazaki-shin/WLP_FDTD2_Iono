#include <iostream>
#include <cmath>
#include "WLP2D.h"

void cal_Ez(double ***Ez,double **Ezp,int q){
  double *f = new double[P+1];
  for(int n = 0; n <= N; n++){
    Fai(s*n*Ez_dt,f);
    for(int i = 0; i < Nx+1; i++){
      for(int j = 0; j < Ny+1; j++){
	Ez[i][j][n] += Ezp[i][j]*f[q];
      }	    
    }
  }
  delete[] f;
}

