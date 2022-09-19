#include <iostream>
#include <cmath>
#include "WLP2D.h"

void cal_Hyp(double **Hyp, double **SUM_Hyp,double **Ezp){

  int ix {(int)(Rx/Dx)-L};
  int jh {(int)(hi/Dy)+L};
  
   /*TE・TMの範囲*/
  for(int i = L; i < ix; i++){
    for(int j = L+1; j <= jh; j++){
      Hyp[i][j] = b/Dx*(Ezp[i+1][j] - Ezp[i][j]) - 2*SUM_Hyp[i][j];
    }
  }

 /*TE・TM_PMLの範囲*/
  double *btx0 =  memory_allocate1d(Nx+1, 0.0); /* σx */
  double *bty0 =  memory_allocate1d(Ny+1, 0.0);/* σy */
  double *btx1 =  memory_allocate1d(Nx, 0.0); /* σx* */
  double *bty1 =  memory_allocate1d(Ny, 0.0); /* σy* */
  
  PML_cond(btx0, btx1, bty0, bty1);
  
  /* (1,1)~(L,jh) */
  for(int i = 0; i < L; i++){
    for(int j = 1; j <= jh; j++){

      Hyp[i][j] = btx1[i]/MU0/Dx*(Ezp[i+1][j] - Ezp[i][j]) - s*btx1[i]*SUM_Hyp[i][j];
    }
  }

 /* (ix,1)~(ie,jh) */

   int ie {(int)(Rx/Dx)-1};
   
   for(int i = ix; i <= ie; i++){
     for(int j = 1; j <= jh; j++){
   
       Hyp[i][j] = btx1[i]/MU0/Dx*(Ezp[i+1][j] - Ezp[i][j]) - s*btx1[i]*SUM_Hyp[i][j];

     }
   }

   /* (L+1,1)~(ix,L) */
   
   for(int i = L; i < ix; i++){
     for(int j = 1; j <= L; j++){
       
      Hyp[i][j] = btx1[i]/MU0/Dx*(Ezp[i+1][j] - Ezp[i][j]) - s*btx1[i]*SUM_Hyp[i][j];
       
     }
   }

    /*Ionoの範囲*/

   int jb {(int)(Ry/Dy)-L};
   
   for(int i = L; i < ix; i++){
     for(int j = jh+1; j < jb; j++){
       Hyp[i][j] = 2/MU0/s/Dx*(Ezp[i+1][j] - Ezp[i][j]) - 2*SUM_Hyp[i][j];
     }
   }

   /*Iono_PMLの範囲*/

    int jp {(int)(Ry/Dy)-1};
    int jq {(int)(Ry/Dy)-L};
    
   /* (1,jh+1) ~ (L,jp)*/
    for(int i = 0; i < L; i++){
      for(int j = jh+1; j <= jp; j++){
	Hyp[i][j] = btx1[i]/MU0/Dx*(Ezp[i+1][j] - Ezp[i][j]) - s*btx1[i]*SUM_Hyp[i][j];
    } 
  }

     /* (Nx-L,jh+1) ~ (ie,jp)*/
    for(int i = ix; i <= ie; i++){
      for(int j = jh+1; j <= jp; j++){
	
	Hyp[i][j] = btx1[i]/MU0/Dx*(Ezp[i+1][j] - Ezp[i][j]) - s*btx1[i]*SUM_Hyp[i][j];
      }
    }

    /* (L+1,jq) ~ (ix,jp)*/
    for(int i = L; i < ix; i++){
      for(int j = jq; j <= jp; j++){
	Hyp[i][j] = btx1[i]/MU0/Dx*(Ezp[i+1][j] - Ezp[i][j]) - s*btx1[i]*SUM_Hyp[i][j];
      }
    }
}
