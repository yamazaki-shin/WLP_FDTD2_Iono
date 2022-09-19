#include <iostream>
#include <cmath>
#include "WLP2D.h"

void cal_Jmxp(double **Jmxp, double **SUM_Jmxp, double **SUM_Jmyp, double **SUM_Jmzp, double **Exp, double **Eyp, double **Ezp){

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
   
   for(int i = L; i < ix; i++){
     for(int j = jh+1; j < ib; j++){
       Jmxp[i][j] = 1/wc*((s1*s1 + bx*bx)*(EPS0*wp*wp*Exp[i][j] - s*SUM_Jmxp[i][j]) + (bx*by - s1*bz)*(EPS0*wp*wp*Eyp[i][j] - s*SUM_Jmyp[i][j]) + (by*s1 + bz*bx)*(EPS0*wp*wp*Ezp[i][j] - s*SUM_Jmzp[i][j]));
     }
   }

}
