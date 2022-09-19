#include <iostream>
#include <cmath>
#include "WLP2D.h"

void cal_Hxp(double **Hxp, double **SUM_Hxp,double **Ezp){

  int ix {(int)(Rx/Dx)-L};
  int jh {(int)(hi/Dy)+L};
  
    /*TE・TMの範囲*/
  for(int i = L+1; i < ix; i++){
    for(int j = L; j < jh; j++){
      Hxp[i][j] = -b/Dy*(Ezp[i][j+1] - Ezp[i][j]) - 2*SUM_Hxp[i][j];
    }
  }


  /*TE・TM_PMLの範囲*/
  double *btx0 =  memory_allocate1d(Nx+1, 0.0); /* σx */
  double *bty0 =  memory_allocate1d(Ny+1, 0.0);/* σy */
  double *btx1 =  memory_allocate1d(Nx, 0.0); /* σx* */
  double *bty1 =  memory_allocate1d(Ny, 0.0); /* σy* */
  
  PML_cond(btx0, btx1, bty0, bty1);
  
  /* (1,1)~(L,jh) */
  for(int i = 1; i <= L; i++){
    for(int j = 1; j < jh; j++){

      Hxp[i][j] = -bty1[j]/MU0/Dy*(Ezp[i][j+1] - Ezp[i][j]) - s*bty1[j]*SUM_Hxp[i][j];
    }
  }

 /* (ix,1)~(ie,jh) */

   int ie {(int)(Rx/Dx)-1};
   
   for(int i = ix; i <= ie; i++){
     for(int j = 1; j < jh; j++){
   
       Hxp[i][j] = -bty1[j]/MU0/Dy*(Ezp[i][j+1] - Ezp[i][j]) - s*bty1[j]*SUM_Hxp[i][j];

     }
   }

   /* (L+1,1)~(ix,L) */
   
   for(int i = L+1; i < ix; i++){
     for(int j = 1; j < L; j++){
       
       Hxp[i][j] = -bty1[j]/MU0/Dy*(Ezp[i][j+1] - Ezp[i][j]) - s*bty1[j]*SUM_Hxp[i][j];
       
     }
   }

   /*Ionoの範囲*/

   int ib {(int)(Ry/Dy)-L};
   
   for(int i = L+1; i < ix; i++){
     for(int j = jh; j < ib; j++){
       Hxp[i][j] = -2/MU0/s/Dy*(Ezp[i][j+1] - Ezp[i][j]) - 2*SUM_Hxp[i][j];
     }
   }

   /*Iono_PMLの範囲*/

    int jp {(int)(Ry/Dy)-1};
    int jq {(int)(Ry/Dy)-L};
    
   /* (1,jh+1) ~ (L,jp)*/
    for(int i = 1; i <= L; i++){
      for(int j = jh; j <= jp; j++){
	Hxp[i][j] = -bty1[j]/MU0/Dy*(Ezp[i][j+1] - Ezp[i][j]) - s*bty1[j]*SUM_Hxp[i][j];
    } 
  }

     /* (Nx-L,jh+1) ~ (ie,jp)*/
    for(int i = ix; i <= ie; i++){
      for(int j = jh; j <= jp; j++){
	
	Hxp[i][j] = -bty1[j]/MU0/Dy*(Ezp[i][j+1] - Ezp[i][j]) - s*bty1[j]*SUM_Hxp[i][j];
      }
    }

    /* (L+1,jq) ~ (ix,jp)*/
    for(int i = L+1; i < ix; i++){
      for(int j = jq; j <= jp; j++){
	Hxp[i][j] = -bty1[j]/MU0/Dy*(Ezp[i][j+1] - Ezp[i][j]) - s*bty1[j]*SUM_Hxp[i][j];
      }
    }
}
