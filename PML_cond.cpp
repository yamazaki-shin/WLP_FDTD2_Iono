#include <iostream>
#include <cmath>
#include "WLP2D.h"


void PML_cond(double *btx0, double *btx1, double *bty0, double *bty1){
  
  constexpr double R {1.0e-6};
  constexpr double M {3.4};
  constexpr double S_MAX {-(M + 1)*sqrt(EPS0/MU0)*std::log(R)/2/L/Dx};

  double Sx0;
  double Sy0;
  double Sx1;
  double Sy1;
    
  for(int i=Nx-L; i <= Nx; i++){
    Sx0 = S_MAX*std::pow((i-Nx+L)*Dx/L/Dx,M);
    btx0[i] = 1/(s/2+Sx0);
  }

  for(int i=Nx-L; i < Nx; i++){
    Sx1 = S_MAX*std::pow((i+0.5-Nx+L)*Dx/L/Dx,M)*MU0/EPS0;
    btx1[i] = 1/(s/2+Sx1);
  }
  
  for(int i=0; i <= L; i++){
    Sx0 = S_MAX*std::pow((L*Dx-i*Dx)/L/Dx,M);
    btx0[i] = 1/(s/2+Sx0);
  }	
  
  for(int i=0; i < L ; i++){
    Sx1 =  S_MAX*std::pow((L*Dx-(i+0.5)*Dx)/L/Dx,M)*MU0/EPS0;
    btx1[i] = 1/(s/2+Sx1);
  }
  
  for(int j=Ny-L; j <= Ny; j++){
    Sy0 = S_MAX*std::pow((j-Ny+L)*Dy/L/Dy,M);
    bty0[j] = 1/(s/2+Sy0);
  }

  for(int j=Ny-L; j < Ny; j++){
    Sy1 = S_MAX*std::pow((j+0.5-Ny+L)*Dy/L/Dy,M)*MU0/EPS0;
    bty1[j] = 1/(s/2+Sy1);
  }
  
  for(int j=0; j <= L; j++){
    Sy0 = S_MAX*std::pow((L*Dy-j*Dy)/L/Dy,M);
    bty0[j] = 1/(s/2+Sy0);
  }	
  
  for(int j=0; j < L ; j++){
    Sy1 =  S_MAX*std::pow((L*Dy-(j+0.5)*Dy)/L/Dy,M)*MU0/EPS0;
    bty1[j] = 1/(s/2+Sy1);
  }

}
