#include <iostream>
#include <cmath>
#include <eigen3/Eigen/Sparse>
#include "WLP2D.h"

void cal_RHS_vector(Eigen::VectorXd &b, double **SUM_Hxp, double **SUM_Hyp, double **SUM_Hzp, double **SUM_Hzxp, double **SUM_Hzyp, double **SUM_Exp, double **SUM_Eyp, double **SUM_Ezp,double **SUM_Fxp, double **SUM_Fyp, double **SUM_Fzxp, double **SUM_Fzyp, double **SUM_Jmxp, double **SUM_Jmyp, double **SUM_Jmzp, double *Jzp, int q, int ***nd){

  int ix {(int)(Rx/Dx)-L};
  int jh {(int)(hi/Dy)+L};
  
  /*TE・TMの範囲*/
  for(int i = L; i < ix; i++){       
    for(int j = L+1; j <= jh; j++){
      
      b( nd[EX][i][j] ) = -2*a/Dy*(SUM_Hzp[i][j] - SUM_Hzp[i][j-1]) - a*EPS0*s*SUM_Exp[i][j]; /*Lx*/
      
    }
  }

  for(int i = L+1; i < ix ; i++){      
    for(int j = L; j < jh; j++){
      
      b( nd[EY][i][j] ) = 2*a/Dx*(SUM_Hzp[i][j] - SUM_Hzp[i-1][j]) - a*EPS0*s*SUM_Eyp[i][j]; /*Ly*/
      
    }
  }

  for(int i = L+1; i < ix ; i++){      
    for(int j = L+1; j <= jh; j++){
      
      b( nd[EZ][i][j] ) = 2*a*((SUM_Hxp[i][j] - SUM_Hxp[i][j-1])/Dy - (SUM_Hyp[i][j] - SUM_Hyp[i-1][j])/Dx) - a*EPS0*s*SUM_Ezp[i][j]; /*Lz*/
      
    }
  }

  int js {(int)(Sy/Dy)};
  int is {(int)(Sx/Dx)};
  
  b( nd[EZ][is][js] ) = b( nd[EZ][is][js] ) - a*Jzp[q];



  
  /*TE・TM_PMLの範囲*/
  double *btx0 =  memory_allocate1d(Nx+1, 0.0); /* σx */
  double *bty0 =  memory_allocate1d(Ny+1, 0.0);/* σy */
  double *btx1 =  memory_allocate1d(Nx, 0.0); /* σx* */
  double *bty1 =  memory_allocate1d(Ny, 0.0); /* σy* */
  
  PML_cond(btx0, btx1, bty0, bty1);

  
  /* (1,1)~(L,jh) */
  
  for(int i = 0; i < L; i++){               
    for(int j = 1; j <= jh; j++){

      b( nd[EX][i][j] ) = gm*((2.0-s*bty0[j])*SUM_Fxp[i][j] - 2*SUM_Exp[i][j] - s*bty0[j]/EPS0/Dy*(btx1[i]*(SUM_Hzxp[i][j] - SUM_Hzxp[i][j-1]) + bty1[j]*SUM_Hzyp[i][j] - bty1[j-1]*SUM_Hzyp[i][j-1])); /*Lx*/
       
    }
  }

  
  for(int i = 1; i <= L; i++){               
    for(int j = 0; j < jh; j++){

      b( nd[EY][i][j] ) = gm*((2.0-s*btx0[i])*SUM_Fyp[i][j] - 2*SUM_Eyp[i][j] + s*btx0[i]/EPS0/Dx*(btx1[i]*SUM_Hzxp[i][j] - btx1[i-1]*SUM_Hzxp[i-1][j] + bty1[j]*(SUM_Hzyp[i][j] - SUM_Hzyp[i-1][j]))); /*Ly*/
       
    }
  }

 for(int i = 1; i <= L; i++){
    for(int j = 1; j <= jh; j++){

      b( nd[EZ][i][j] ) = gm*((2.0-s*btx0[i])*SUM_Fzxp[i][j] + (2.0-s*bty0[j])*SUM_Fzyp[i][j] - 2*SUM_Ezp[i][j] - s/EPS0*(-bty0[j]/Dy*(bty1[j]*SUM_Hxp[i][j] - bty1[j-1]*SUM_Hxp[i][j-1]) + btx0[i]/Dx*(btx1[i]*SUM_Hyp[i][j] - btx1[i-1]*SUM_Hyp[i-1][j]))); /*Lz*/
       
    }
  }





 
 /* (ix,1)~(ie,jh) */

 int ie {(int)(Rx/Dx)-1};
  
 for(int i = ix; i <= ie; i++){               
    for(int j = 1; j <= jh; j++){

      b( nd[EX][i][j] ) = gm*((2.0-s*bty0[j])*SUM_Fxp[i][j] - 2.0*SUM_Exp[i][j] - s*bty0[j]/EPS0/Dy*(btx1[i]*(SUM_Hzxp[i][j] - SUM_Hzxp[i][j-1]) + bty1[j]*SUM_Hzyp[i][j] - bty1[j-1]*SUM_Hzyp[i][j-1])); /*Lx*/
       
    }
  }

  
  for(int i = ix; i <= ie; i++){               
    for(int j = 0; j < jh; j++){

      b( nd[EY][i][j] ) = gm*((2.0-s*btx0[i])*SUM_Fyp[i][j] - 2.0*SUM_Eyp[i][j] + s*btx0[i]/EPS0/Dx*(btx1[i]*SUM_Hzxp[i][j] - btx1[i-1]*SUM_Hzxp[i-1][j] + bty1[j]*(SUM_Hzyp[i][j] - SUM_Hzyp[i-1][j]))); /*Ly*/
       
    }
  }

 for(int i = ix; i <= ie; i++){
    for(int j = 1; j <= jh; j++){

        b( nd[EZ][i][j] ) = gm*((2.0-s*btx0[i])*SUM_Fzxp[i][j] + (2.0-s*bty0[j])*SUM_Fzyp[i][j] - 2*SUM_Ezp[i][j] - s/EPS0*(-bty0[j]/Dy*(bty1[j]*SUM_Hxp[i][j] - bty1[j-1]*SUM_Hxp[i][j-1]) + btx0[i]/Dx*(btx1[i]*SUM_Hyp[i][j] - btx1[i-1]*SUM_Hyp[i-1][j]))); /*Lz*/
       
    }
  }






 
/* (L+1,1)~(ix,L) */
 
  for(int i = L; i < ix; i++){               
    for(int j = 1; j <= L; j++){

      b( nd[EX][i][j] ) = gm*((2.0-s*bty0[j])*SUM_Fxp[i][j] - 2.0*SUM_Exp[i][j] - s*bty0[j]/EPS0/Dy*(btx1[i]*(SUM_Hzxp[i][j] - SUM_Hzxp[i][j-1]) + bty1[j]*SUM_Hzyp[i][j] - bty1[j-1]*SUM_Hzyp[i][j-1])); /*Lx*/
       
    }
  }

  
  for(int i = L+1; i < ix; i++){               
    for(int j = 0; j < L; j++){

      b( nd[EY][i][j] ) = gm*((2.0-s*btx0[i])*SUM_Fyp[i][j] - 2.0*SUM_Eyp[i][j] + s*btx0[i]/EPS0/Dx*(btx1[i]*SUM_Hzxp[i][j] - btx1[i-1]*SUM_Hzxp[i-1][j] + bty1[j]*(SUM_Hzyp[i][j] - SUM_Hzyp[i-1][j]))); /*Ly*/
       
    }
  }

 for(int i = L+1; i < ix; i++){
    for(int j = 1; j <= L; j++){

       b( nd[EZ][i][j] ) = gm*((2.0-s*btx0[i])*SUM_Fzxp[i][j] + (2.0-s*bty0[j])*SUM_Fzyp[i][j] - 2*SUM_Ezp[i][j] - s/EPS0*(-bty0[j]/Dy*(bty1[j]*SUM_Hxp[i][j] - bty1[j-1]*SUM_Hxp[i][j-1]) + btx0[i]/Dx*(btx1[i]*SUM_Hyp[i][j] - btx1[i-1]*SUM_Hyp[i-1][j]))); /*Lz*/
    }
  }




 
  
 /*Ionoの範囲*/
 
 double v {0.0};
 cal_Nyu(v,z);
 double s1 {(s/2 + v)/wc};
 double bx {1/sqrt(2)};
 double by {-1/sqrt(2)};
 double bz {0.0};
 
 int jb {(int)(Ry/Dy)-L};
 
 for(int i = L; i < ix; i++){
   for(int j = jh+1; j < jb; j++){
     b( nd[EX][i][j] ) = -4/EPS0/s/Dy*(SUM_Hzp[i][j] -SUM_Hzp[i][j-1]) - 2*SUM_Exp[i][j] - 2/EPS0/wc*((s1*s1+bx*bx)*SUM_Jmxp[i][j] + (bx*by-s1*bz)*SUM_Jmyp[i][j] + (by*s1+bz*bx)*SUM_Jmzp[i][j]); /* Lx */
   }
 }
 
 for(int i = L+1; i < ix; i++){
   for(int j = jh; j < jb; j++){
     b( nd[EY][i][j] ) = 4/EPS0/s/Dx*(SUM_Hzp[i][j] - SUM_Hzp[i-1][j]) - 2*SUM_Eyp[i][j] - 2/EPS0/wc*((bz*s1+bx*by)*SUM_Jmxp[i][j] + (s1*s1+by*by)*SUM_Jmyp[i][j] + (by*bz-s1*bx)*SUM_Jmzp[i][j]); /* Ly */
   }
 }
 
 
 for(int i = L+1; i < ix; i++){
   for(int j = jh + 1; j < jb; j++){
     b( nd[EZ][i][j] ) = 4/EPS0/s*((SUM_Hxp[i][j] - SUM_Hxp[i][j-1])/Dy - (SUM_Hyp[i][j] - SUM_Hyp[i-1][j])/Dx) - 2*SUM_Ezp[i][j] - 2/EPS0/wc*((bz*bx-by*s1)*SUM_Jmxp[i][j] + (bx*s1+by*bz)*SUM_Jmyp[i][j] + (s1*s1-bz*bz)*SUM_Jmzp[i][j]); /* Lz */
   }
 }
 






 
 /*Iono_PMLの範囲*/
 
 int jp {(int)(Ry/Dy)-1};
 int jq {(int)(Ry/Dy)-L};
 
 /* (1,jh+1) ~ (L,jp)*/
 for(int i = 0; i < L; i++){
   for(int j = jh + 1; j <= jp; j++){
     b( nd[EX][i][j] ) = (2.0-s*bty0[j])*SUM_Fxp[i][j] - s*bty0[j]/EPS0/Dy*(btx1[i]*(SUM_Hzxp[i][j] - SUM_Hzxp[i][j-1]) + bty1[j]*SUM_Hzyp[i][j] - bty1[j-1]*SUM_Hzyp[i][j-1]) - 2*SUM_Exp[i][j] + 2/EPS0/wc*((s1*s1+bx*bx)*SUM_Jmxp[i][j] + (bx*by-s1*bz)*SUM_Jmyp[i][j] + (by*s1+bz*bx)*SUM_Jmzp[i][j]); /*Lx*/
     
   }
 }
 
 
 for(int i = 1; i <= L; i++){
   for(int j = jh; j <= jp; j++){
     b( nd[EY][i][j] ) = (2.0-s*btx0[i])*SUM_Fyp[i][j] + s*btx0[i]/EPS0/Dx*(btx1[i]*SUM_Hzxp[i][j] - btx1[i-1]*SUM_Hzxp[i-1][j] + bty1[j]*(SUM_Hzyp[i][j] - SUM_Hzyp[i-1][j])) - 2*SUM_Eyp[i][j] + 2/EPS0/wc*((bz*s1+bx*by)*SUM_Jmxp[i][j] + (s1*s1+by*by)*SUM_Jmyp[i][j] + (by*bz-s1*bx)*SUM_Jmzp[i][j]); /*Ly*/
   }
 }
 
 for(int i = 1; i <= L; i++){
   for(int j = jh + 1; j <= jp; j++){
     b( nd[EZ][i][j] ) = (2.0-s*btx0[i])*SUM_Fzxp[i][j] + (2.0-s*bty0[j])*SUM_Fzyp[i][j] - s/EPS0*(-bty0[j]/Dy*(bty1[j]*SUM_Hxp[i][j] - bty1[j-1]*SUM_Hxp[i][j-1]) + btx0[i]/Dx*(btx1[i]*SUM_Hyp[i][j] - btx1[i-1]*SUM_Hyp[i-1][j])) - 2*SUM_Ezp[i][j] + 2/EPS0/wc*((bz*bx-by*s1)*SUM_Jmxp[i][j] + (bx*s1+by*bz)*SUM_Jmyp[i][j] + (s1*s1-bz*bz)*SUM_Jmzp[i][j]); /*Lz*/
    }
  }






 
   /* (Nx-L,jh+1) ~ (ie,jp)*/
   for(int i = ix; i <= ie; i++){
    for(int j = jh + 1; j <= jp; j++){
      b( nd[EX][i][j] ) = (2.0-s*bty0[j])*SUM_Fxp[i][j] - s*bty0[j]/EPS0/Dy*(btx1[i]*(SUM_Hzxp[i][j] - SUM_Hzxp[i][j-1]) + bty1[j]*SUM_Hzyp[i][j] - bty1[j-1]*SUM_Hzyp[i][j-1]) - 2*SUM_Exp[i][j] + 2/EPS0/wc*((s1*s1+bx*bx)*SUM_Jmxp[i][j] + (bx*by-s1*bz)*SUM_Jmyp[i][j] + (by*s1+bz*bx)*SUM_Jmzp[i][j]); /*Lx*/
      
    }
  }
  
  
  for(int i = ix; i <= ie; i++){
    for(int j = jh; j <= jp; j++){
      b( nd[EY][i][j] ) = (2.0-s*btx0[i])*SUM_Fyp[i][j] + s*btx0[i]/EPS0/Dx*(btx1[i]*SUM_Hzxp[i][j] - btx1[i-1]*SUM_Hzxp[i-1][j] + bty1[j]*(SUM_Hzyp[i][j] - SUM_Hzyp[i-1][j])) - 2*SUM_Eyp[i][j] + 2/EPS0/wc*((bz*s1+bx*by)*SUM_Jmxp[i][j] + (s1*s1+by*by)*SUM_Jmyp[i][j] + (by*bz-s1*bx)*SUM_Jmzp[i][j]); /*Ly*/
    }
  }
  
  for(int i = ix; i <= ie; i++){
    for(int j = jh + 1; j <= jp; j++){
   b( nd[EZ][i][j] ) = (2.0-s*btx0[i])*SUM_Fzxp[i][j] + (2.0-s*bty0[j])*SUM_Fzyp[i][j] - s/EPS0*(-bty0[j]/Dy*(bty1[j]*SUM_Hxp[i][j] - bty1[j-1]*SUM_Hxp[i][j-1]) + btx0[i]/Dx*(btx1[i]*SUM_Hyp[i][j] - btx1[i-1]*SUM_Hyp[i-1][j])) - 2*SUM_Ezp[i][j] + 2/EPS0/wc*((bz*bx-by*s1)*SUM_Jmxp[i][j] + (bx*s1+by*bz)*SUM_Jmyp[i][j] + (s1*s1-bz*bz)*SUM_Jmzp[i][j]); /*Lz*/
    }
  }



  


  
 /* (L+1,jq) ~ (ix,jp)*/
  for(int i = L; i < ix; i++){
    for(int j = jq; j <= jp; j++){
      b( nd[EX][i][j] ) = (2.0-s*bty0[j])*SUM_Fxp[i][j] - s*bty0[j]/EPS0/Dy*(btx1[i]*(SUM_Hzxp[i][j] - SUM_Hzxp[i][j-1]) + bty1[j]*SUM_Hzyp[i][j] - bty1[j-1]*SUM_Hzyp[i][j-1]) - 2*SUM_Exp[i][j] + 2/EPS0/wc*((s1*s1+bx*bx)*SUM_Jmxp[i][j] + (bx*by-s1*bz)*SUM_Jmyp[i][j] + (by*s1+bz*bx)*SUM_Jmzp[i][j]); /*Lx*/
      
    }
  }
  
  
  for(int i = L+1; i < ix; i++){
    for(int j = jq; j <= jp; j++){
        b( nd[EY][i][j] ) = (2.0-s*btx0[i])*SUM_Fyp[i][j] + s*btx0[i]/EPS0/Dx*(btx1[i]*SUM_Hzxp[i][j] - btx1[i-1]*SUM_Hzxp[i-1][j] + bty1[j]*(SUM_Hzyp[i][j] - SUM_Hzyp[i-1][j])) - 2*SUM_Eyp[i][j] + 2/EPS0/wc*((bz*s1+bx*by)*SUM_Jmxp[i][j] + (s1*s1+by*by)*SUM_Jmyp[i][j] + (by*bz-s1*bx)*SUM_Jmzp[i][j]); /*Ly*/
    }
  }
  
  for(int i = L+1; i <= ix; i++){
    for(int j = jq; j <= jp; j++){
      b( nd[EZ][i][j] ) = (2.0-s*btx0[i])*SUM_Fzxp[i][j] + (2.0-s*bty0[j])*SUM_Fzyp[i][j] - s/EPS0*(-bty0[j]/Dy*(bty1[j]*SUM_Hxp[i][j] - bty1[j-1]*SUM_Hxp[i][j-1]) + btx0[i]/Dx*(btx1[i]*SUM_Hyp[i][j] - btx1[i-1]*SUM_Hyp[i-1][j])) - 2*SUM_Ezp[i][j] + 2/EPS0/wc*((bz*bx-by*s1)*SUM_Jmxp[i][j] + (bx*s1+by*bz)*SUM_Jmyp[i][j] + (s1*s1-bz*bz)*SUM_Jmzp[i][j]); /*Lz*/
    }
  }
  


  
  delete1d(btx0);
  delete1d(btx1);
  delete1d(bty0);
  delete1d(bty1);
  
}
