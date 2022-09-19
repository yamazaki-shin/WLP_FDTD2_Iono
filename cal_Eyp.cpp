#include <iostream>
#include <cmath>
#include <eigen3/Eigen/Sparse>
#include "WLP2D.h"

void cal_Eyp(double **Eyp, int ***nd, Eigen::VectorXd &b, Eigen::SparseLU < Eigen::SparseMatrix <double>, Eigen::COLAMDOrdering <int> > &S){

  Eigen::VectorXd x = S.solve(b);
  
  for(int i = 0; i <= Nx; i++){
    for(int j = 0; j < Ny; j++){
      Eyp[i][j] = x( nd[EY][i][j] );
    }
  }
}
