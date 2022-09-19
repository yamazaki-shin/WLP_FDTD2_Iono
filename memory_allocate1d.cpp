#include <iostream>
#include <cmath>
#include "WLP2D.h"

double *memory_allocate1d(int M, double initial_value){
  double *v = new double [M];
  for(int n = 0; n < M; n++){
    v[n] = initial_value;
  }
  return v;
}
