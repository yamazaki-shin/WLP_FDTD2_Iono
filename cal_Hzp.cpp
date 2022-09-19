#include <iostream>
#include <cmath>
#include "WLP2D.h"

void cal_Hzp(double **Hzp,double **Hzxp, double **Hzyp, double **SUM_Hzp,double **Exp, double **Eyp){

  int ix {(int)(Rx/Dx)-L};
  int jh {(int)(hi/Dy)+L};
  
  /*TE・TMの範囲*/
  for(int i = L+1; i <= ix; i++){
    for(int j = L+1; j <= jh; j++){
      Hzp[i][j] = -b*(-(Exp[i][j+1] - Exp[i][j])/Dy + (Eyp[i+1][j] - Eyp[i][j])/Dx) - 2*SUM_Hzp[i][j];
    }
  }

/*TE・TM_PMLの範囲*/
   /* (1,1)~(L,jh) */
  for(int i = 1; i <= L; i++){
    for(int j = 1; j <= jh; j++){
      Hzp[i][j] = Hzxp[i][j] + Hzyp[i][j];
    }
  }

  /* (Nx-L,1)~(Nx-1,jh) */

   int ie {(int)(Rx/Dx)-1};
   
   for(int i = ix; i <= ie; i++){
     for(int j = 1; j <= jh; j++){
     Hzp[i][j] = Hzxp[i][j] + Hzyp[i][j];
    }
  }

  /* (L+1,1)~(Nx-L-1,L) */
   
   for(int i = L+1; i < ix; i++){
     for(int j = 1; j <= L; j++){
      Hzp[i][j] = Hzxp[i][j] + Hzyp[i][j];
    }
  }   
  
 /*Ionoの範囲*/

  int ib {(int)(Ry/Dy)-L-1};

  for(int i = L+1; i <= ix; i++){
    for(int j = jh+1; j <= ib; j++){
      Hzp[i][j] = -2/MU0/s*((Eyp[i+1][j] - Eyp[i][j])/Dx - (Exp[i][j+1] - Exp[i][j])/Dy) - 2*SUM_Hzp[i][j];
    }
  }

   /*Iono_PMLの範囲*/

    int jp {(int)(Ry/Dy)-1};
    int jq {(int)(Ry/Dy)-L};
    
   /* (1,jh+1) ~ (L,jp)*/
    for(int i = 1; i < L; i++){
      for(int j = jh+1; j <= jp; j++){
	 Hzp[i][j] = Hzxp[i][j] + Hzyp[i][j];
    }
  }

     /* (Nx-L,jh+1) ~ (ie,jp)*/
    for(int i = ix; i <= ie; i++){
      for(int j = jh+1; j <= jp; j++){
	Hzp[i][j] = Hzxp[i][j] + Hzyp[i][j];
      }
    }	

     /* (L+1,jq) ~ (ix,jp)*/
    for(int i = L+1; i <= ix; i++){
      for(int j = jq; j <= jp; j++){
    	Hzp[i][j] = Hzxp[i][j] + Hzyp[i][j];
      }
    }	
}
