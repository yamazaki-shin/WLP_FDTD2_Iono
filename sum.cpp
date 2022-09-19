#include <iostream>
#include <cmath>
#include "WLP2D.h"

void sum(double **SUM, double **f, int n, int m){
  for(int i = 0; i < n; i++){
    for(int j = 0; j < m; j++){
      SUM[i][j] += f[i][j];
    }
  }
}
