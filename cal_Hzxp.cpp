#include <iostream>
#include <cmath>
#include "WLP2D.h"

void cal_Hzxp(double **Hzxp, double **SUM_Hzxp, double **Eyp){

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
    for(int j = 1; j <= jh; j++){

      Hzxp[i][j] = -btx1[i]/MU0/Dx*(Eyp[i+1][j] - Eyp[i][j]) - s*btx1[i]*SUM_Hzxp[i][j];
    }
  }

 /* (Nx-L,1)~(Nx-1,jh) */

   int ie {(int)(Rx/Dx)-1};
   
   for(int i = ix; i <= ie; i++){
     for(int j = 1; j <= jh; j++){
   
         Hzxp[i][j] = -btx1[i]/MU0/Dx*(Eyp[i+1][j] - Eyp[i][j]) - s*btx1[i]*SUM_Hzxp[i][j];

     }
   }

   /* (L+1,1)~(Nx-L-1,L) */
   
   for(int i = L+1; i < ix; i++){
     for(int j = 1; j <= L; j++){
       
          Hzxp[i][j] = -btx1[i]/MU0/Dx*(Eyp[i+1][j] - Eyp[i][j]) - s*btx1[i]*SUM_Hzxp[i][j];
       
     }
   }


   
   /*Iono_PMLの範囲*/

    int jp {(int)(Ry/Dy)-1};
    int jq {(int)(Ry/Dy)-L};
    
   /* (1,jh+1) ~ (L,jp)*/
    for(int i = 1; i < L; i++){
      for(int j = jh+1; j <= jp; j++){
	Hzxp[i][j] = -btx1[i]/MU0/Dx*(Eyp[i+1][j] - Eyp[i][j]) - s*btx1[i]*SUM_Hzxp[i][j];
    } 
  }

     /* (Nx-L,jh+1) ~ (ie,jp)*/
    for(int i = ix; i <= ie; i++){
      for(int j = jh+1; j <= jp; j++){
	Hzxp[i][j] = -btx1[i]/MU0/Dx*(Eyp[i+1][j] - Eyp[i][j]) - s*btx1[i]*SUM_Hzxp[i][j];	
      }
    }

    /* (L+1,jq) ~ (ix,jp)*/
    for(int i = L+1; i <= ix; i++){
      for(int j = jq; j <= jp; j++){
	Hzxp[i][j] = -btx1[i]/MU0/Dx*(Eyp[i+1][j] - Eyp[i][j]) - s*btx1[i]*SUM_Hzxp[i][j];
      }
    }
   
}
