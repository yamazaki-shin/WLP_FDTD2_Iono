#include <iostream>
#include "WLP2D.h"

void cal_nd(int ***nd){
  int nd_iter = 0;
  
  for(int i = 0; i <= Nx; i++){
    for(int j = 0; j <= Ny; j++){
      
    double  Ex_x = (i + 0.5)*Dx;
    double  Ex_y = j*Dy;
    double  Ey_x = i*Dx;
    double  Ey_y = (j + 0.5)*Dy;
    double  Ez_x = i*Dx;
    double  Ez_y = j*Dy;
   
      
      if((Ex_x <= Rx) && (Ex_y <= Ry)){
	nd[EX][i][j] = nd_iter;
	nd_iter++;
      }
      
      if((Ey_x <= Rx) && (Ey_y <= Ry)){
	nd[EY][i][j] = nd_iter;
	nd_iter++;
      }
      
      if((Ez_x <= Rx) && (Ez_y <= Ry)){
	nd[EZ][i][j] = nd_iter;
	nd_iter++;
      }
     
    }
  }
}
