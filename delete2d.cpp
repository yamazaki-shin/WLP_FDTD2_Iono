#include <iostream>
#include <cmath>
#include "WLP2D.h"

void delete2d(double **v, int M){
  for(int m = 0; m < M; m++){
    delete [] v[m];
  }
  delete [] v;
}
