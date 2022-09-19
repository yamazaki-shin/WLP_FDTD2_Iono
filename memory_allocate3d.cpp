#include <iostream>
#include <cmath>
#include "WLP2D.h"

double ***memory_allocate3d(int M, int N, int L,  double initial_value){
  double ***v = new double** [M];
  for(int m = 0; m < M; m++){
    v[m] = new double* [N];
    for(int n = 0; n < N; n++){
      v[m][n] = new double [L];
    }
  }
    for(int m = 0; m < M; m++){
      for(int n = 0; n < N; n++){
	for(int l = 0; l < L; l++){     
	  v[m][n][l] = initial_value;
	}
      }
    }
    return v;
}
