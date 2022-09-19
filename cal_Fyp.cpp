#include <iostream>
#include <cmath>
#include "WLP2D.h"

void cal_Fyp(double **Fyp, double **SUM_Fyp,double **Hzp){

  int ix {(int)(Rx/Dx)-L};
  int jh {(int)(hi/Dy)+L};
  
 
 /*TE・TM_PMLの範囲*/
  double *btx0 =  memory_allocate1d(Nx+1, 0.0); /* σx */
  double *bty0 =  memory_allocate1d(Ny+1, 0.0);/* σy */
  double *btx1 =  memory_allocate1d(Nx, 0.0); /* σx* */
  double *bty1 =  memory_allocate1d(Ny, 0.0); /* σy* */
  
  PML_cond(btx0, btx1, bty0, bty1);
  
  /* (1,1)~(L,jh) */
  for(int i = 1; i <= L; i++){
    for(int j = 1; j < jh; j++){

   	Fyp[i][j] = -btx0[i]/EPS0/Dx*(Hzp[i][j] - Hzp[i-1][j]) - s*btx0[j]*SUM_Fyp[i][j];
    }
  }

 /* (ix,1)~(ie,jh) */

   int ie {(int)(Rx/Dx)-1};
   
   for(int i = ix; i <= ie; i++){
     for(int j = 1; j < jh; j++){
   
    	Fyp[i][j] = -btx0[i]/EPS0/Dx*(Hzp[i][j] - Hzp[i-1][j]) - s*btx0[j]*SUM_Fyp[i][j];

     }
   }

   /* (L+1,1)~(ix,L) */
   
   for(int i = L+1; i < ix; i++){
     for(int j = 1; j < L; j++){
       
  	Fyp[i][j] = -btx0[i]/EPS0/Dx*(Hzp[i][j] - Hzp[i-1][j]) - s*btx0[j]*SUM_Fyp[i][j];
       
     }
   }

  /*Iono_PMLの範囲*/

    int jp {(int)(Ry/Dy)-1};
    int jq {(int)(Ry/Dy)-L};
    
   /* (1,jh+1) ~ (L,jp)*/
    for(int i = 1; i <= L; i++){
      for(int j = jh; j <= jp; j++){
	Fyp[i][j] = -btx0[i]/EPS0/Dx*(Hzp[i][j] - Hzp[i-1][j]) - s*btx0[j]*SUM_Fyp[i][j];
    } 
  }

     /* (Nx-L,jh+1) ~ (ie,jp)*/
    for(int i = ix; i <= ie; i++){
      for(int j = jh; j <= jp; j++){
	
	Fyp[i][j] = -btx0[i]/EPS0/Dx*(Hzp[i][j] - Hzp[i-1][j]) - s*btx0[j]*SUM_Fyp[i][j];
      }
    }

    /* (L+1,jq) ~ (ix,jp)*/
    for(int i = L+1; i < ix; i++){
      for(int j = jq; j <= jp; j++){
	Fyp[i][j] = -btx0[i]/EPS0/Dx*(Hzp[i][j] - Hzp[i-1][j]) - s*btx0[j]*SUM_Fyp[i][j];
      }
    }
}
