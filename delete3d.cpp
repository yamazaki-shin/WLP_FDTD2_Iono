#include <iostream>
#include <cmath>
#include "WLP2D.h"

void delete3d(double ***v, int M, int N){
  for(int m = 0; m < M; m++){
    for(int n = 0; n < N; n++){
      delete [] v[m][n];
    }
    delete [] v[m];
  }
  delete [] v;
}
