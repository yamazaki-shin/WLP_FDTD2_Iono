#include <iostream>
#include <cmath>
#include <eigen3/Eigen/Sparse>
#include "WLP2D.h"

void cal_Ezp(double **Ezp, int ***nd, Eigen::VectorXd &b, Eigen::SparseLU < Eigen::SparseMatrix <double>, Eigen::COLAMDOrdering <int> > &S){

  Eigen::VectorXd x = S.solve(b);
  
  for(int i = 0; i <= Nx; i++){
    for(int j = 0; j <= Ny; j++){
      Ezp[i][j] = x( nd[EZ][i][j] );
    }
  }
}
