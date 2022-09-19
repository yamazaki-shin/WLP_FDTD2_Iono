#include <iostream>
#include <cmath>
#include "WLP2D.h"

void cal_P(double **P){
  double v {0};
  cal_Nyu(v,z);
  double s1 {(s/2 + v)/wc};
  double bx {1/sqrt(2)};
  double by {-1/sqrt(2)};
  double bz {0};
  
  P[0][0] = 2/s*wp*wp/wc*(s1*s1 + bx*bx);
  P[0][1] = 2/s*wp*wp/wc*(bx*by - s1*bz);
  P[0][2] = 2/s*wp*wp/wc*(by*s1 + bz*bx);
  P[1][0] = 2/s*wp*wp/wc*(bz*s1 + bx*by);
  P[1][1] = 2/s*wp*wp/wc*(s1*s1 + by*by);
  P[1][2] = 2/s*wp*wp/wc*(by*bz - bx*s1);
  P[2][0] = 2/s*wp*wp/wc*(bz*bx - by*s1);
  P[2][1] = 2/s*wp*wp/wc*(bx*s1 + by*bz);
  P[2][2] = 2/s*wp*wp/wc*(s1*s1 + bz*bz);
  
    }
