#include <iostream>
#include <cmath>
#include "WLP2D.h"

void cal_Jmyp(double **Jmyp, double **SUM_Jmxp, double **SUM_Jmyp, double **SUM_Jmzp, double **Exp, double **Eyp, double **Ezp){

  double v {0.0};
  cal_Nyu(v,z);
  double s1 {(s/2 + v)/wc};
  double bx {1/sqrt(2)};
  double by {-1/sqrt(2)};
  double bz {0};
  
  int ix {(int)(Rx/Dx)-L};
  int jh {(int)(hi/Dy)+L};
  
   /*Ionoの範囲*/

   int ib {(int)(Ry/Dy)-L-1};
   
   for(int i = L+1; i < ix; i++){
     for(int j = jh; j < ib; j++){
       Jmyp[i][j] = 1/wc*((bz*s1 + bx*by)*(EPS0*wp*wp*Exp[i][j] - s*SUM_Jmxp[i][j]) + (s1*s1 + by*by)*(EPS0*wp*wp*Eyp[i][j] - s*SUM_Jmyp[i][j]) + (by*bz - bx*s1)*(EPS0*wp*wp*Ezp[i][j] - s*SUM_Jmzp[i][j]));
     }
   }

}
