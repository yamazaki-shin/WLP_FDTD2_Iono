#include <iostream>
#include <cmath>
#include "WLP2D.h"

int ***memory_allocate3d_int(int M, int N, int L,  int initial_value){
  int ***v = new int** [M];
  for(int m = 0; m < M; m++){
    v[m] = new int* [N];
    for(int n = 0; n < N; n++){
      v[m][n] = new int [L];
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
