#include <iostream>
#include <cmath>
#include "WLP2D.h"

double **memory_allocate2d(int M, int N, double initial_value){
  double **v = new double* [M];
  for(int m = 0; m < M; m++){
    v[m] = new double [N];
    for(int n = 0; n < N; n++){
      v[m][n] = initial_value;
    }
  }
  return v;
}
