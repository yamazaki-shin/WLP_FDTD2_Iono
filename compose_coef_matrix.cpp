#include <iostream>
#include <cmath>
#include <vector>
#include "WLP2D.h"

void compose_coef_matrix(Eigen::SparseLU < Eigen::SparseMatrix <double>, Eigen::COLAMDOrdering <int> > &S, int ***nd, double **P){
  
  std::vector <Eigen::Triplet <double> > list;

  int ix {(int)(Rx/Dx)-L};
  int jh {(int)(hi/Dy)+L};
 
   
  /*TE・TMの範囲*/
  for(int i = L; i < ix ; i++){     
    for(int j = L+1; j <= jh; j++){
      
      /*Lx: Ex^p(i,j), Ex^p(i,j-1), Ex^p(i,j+1) */
      list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EX][i][j] , 1+2*a*b/Dy/Dy));
      list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EX][i][j+1], -a*b/Dy/Dy));
      list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EX][i][j-1], -a*b/Dy/Dy));

      /*Lx: Ey^p(i,j), Ey^p(i,j-1), Ey^p(i+1,j), Ey^p(i+1,j-1) */
      list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EY][i][j], -a*b/Dx/Dy) );
      list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EY][i][j-1], a*b/Dx/Dy) );
      list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EY][i+1][j], a*b/Dx/Dy) );
      list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EY][i+1][j-1], -a*b/Dx/Dy) );
      
    }
  }

  for(int i = L+1; i < ix ; i++){    
    for(int j = L; j < jh; j++){
      /*Ly: Ex^p(i,j), Ex^p(i,j+1), Ex^p(i-1,j+1), Ex^p(i-1,j) */
      list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EX][i][j], -a*b/Dx/Dy));
      list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EX][i][j+1], a*b/Dx/Dy));
      list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EX][i-1][j+1], -a*b/Dx/Dy));
      list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EX][i-1][j], a*b/Dx/Dy));

      /*Ly: Ey^p(i,j), Ey^p(i+1,j), Ey^p(i-1,j) */
      list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EY][i][j], 1+2*a*b/Dx/Dx) );
      list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EY][i+1][j], -a*b/Dx/Dx) );
      list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EY][i-1][j], -a*b/Dx/Dx) );
    }
  }
 
  for(int i = L+1; i < ix ; i++){    
    for(int j = L+1; j <= jh; j++){
      
      /*Lz: Ez^p(i,j), Ez^p(i,j+1), Ez^p(i,j-1), Ez^p(i+1,j), Ez^p(i-1,j) */
      list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i][j], 1+2*a*b/Dx/Dx+2*a*b/Dy/Dy) );
      list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i][j+1], -a*b/Dy/Dy) );
      list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i][j-1], -a*b/Dy/Dy) );
      list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i+1][j], -a*b/Dx/Dx) );
      list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i-1][j], -a*b/Dx/Dx) );
      
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
       
       /*Lx: Ex^p(i,j), Ex^p(i,j+1), Ex^p(i,j-1) */
       list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EX][i][j] , 1+C*C*gm*bty0[j]/Dy/Dy*(bty1[j]+bty1[j-1])));
       list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EX][i][j+1], -C*C*gm*bty0[j]*bty1[j]/Dy/Dy));
       list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EX][i][j-1], -C*C*gm*bty0[j]*bty1[j-1]/Dy/Dy));

        /*Lx: Ey^p(i,j), Ey^p(i,j-1), Ey^p(i+1,j), Ey^p(i+1,j-1) */
       list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EY][i][j], -C*C*gm*bty0[j]*btx1[i]/Dx/Dy) );
       list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EY][i][j-1],  C*C*gm*bty0[j]*btx1[i]/Dx/Dy) );
       list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EY][i+1][j],  C*C*gm*bty0[j]*btx1[i]/Dx/Dy) );
       list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EY][i+1][j-1],  -C*C*gm*bty0[j]*btx1[i]/Dx/Dy) );
     }
   }
   
   for(int i = 1; i <= L; i++){
     for(int j = 0; j < jh; j++){
       
       /*Ly: Ex^p(i,j), Ex^p(i,j+1), Ex^p(i-1,j+1), Ex^p(i-1,j) */
       list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EX][i][j], -C*C*gm*bty1[j]*btx0[i]/Dx/Dy));
       list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EX][i][j+1], C*C*gm*bty1[j]*btx0[i]/Dx/Dy));
       list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EX][i-1][j+1], -C*C*gm*bty1[j]*btx0[i]/Dx/Dy)); 
       list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EX][i-1][j], C*C*gm*bty1[j]*btx0[i]/Dx/Dy));

       /*Ly: Ey^p(i,j), Ey^p(i+1,j), Ey^p(i-1,j) */
       list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EY][i][j],1+C*C*gm*btx0[i]/Dx/Dx*(btx1[i]+btx1[i-1]) ));
       list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EY][i+1][j], -C*C*gm*btx0[i]*btx1[i]/Dx/Dx) );
       list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EY][i-1][j],  -C*C*gm*btx0[i]*btx1[i-1]/Dx/Dx) );
     }
   }
   
   
   for(int i = 1; i <= L; i++){
     for(int j = 1; j <= jh; j++){
       
       /*Lz: Ez^p(i,j), Ez^p(i,j+1), Ez^p(i,j-1), Ez^p(i+1,j), Ez^p(i-1,j) */
       list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i][j], 1 + C*C*gm*bty0[j]/Dy/Dy*(bty1[j]+bty1[j-1]) + C*C*gm*btx0[i]/Dx/Dx*(btx1[i]+btx1[i-1]) ) );
       list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i][j+1],  -C*C*gm*bty0[j]*bty1[j]/Dy/Dy));
       list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i][j-1],  -C*C*gm*bty0[j]*bty1[j-1]/Dy/Dy));
       list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i+1][j],  -C*C*gm*btx1[i]*btx0[i]/Dx/Dx));
       list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i-1][j],  -C*C*gm*btx1[i-1]*btx0[i]/Dx/Dx));
       
     }
   }
   


   
   
   /* (Nx-L,1)~(Nx-1,jh) */
   
   int ie {(int)(Rx/Dx)-1};
   
   for(int i = ix; i <= ie; i++){
     for(int j = 1; j <= jh; j++){
       
       /*Lx: Ex^p(i,j), Ex^p(i,j+1), Ex^p(i,j-1) */
       list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EX][i][j] , 1+C*C*gm*bty0[j]/Dy/Dy*(bty1[j]+bty1[j-1])));
       list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EX][i][j+1], -C*C*gm*bty0[j]*bty1[j]/Dy/Dy));
       list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EX][i][j-1], -C*C*gm*bty0[j]*bty1[j-1]/Dy/Dy));

       /*Lx: Ey^p(i,j), Ey^p(i,j-1), Ey^p(i+1,j), Ey^p(i+1,j-1) */
       list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EY][i][j], -C*C*gm*bty0[j]*btx1[i]/Dx/Dy) );
       list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EY][i][j-1],  C*C*gm*bty0[j]*btx1[i]/Dx/Dy) );
       list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EY][i+1][j],  C*C*gm*bty0[j]*btx1[i]/Dx/Dy) );
       list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EY][i+1][j-1],  -C*C*gm*bty0[j]*btx1[i]/Dx/Dy) );
     }
   }
   
   for(int i = ix; i <= ie; i++){
     for(int j = 0; j < jh; j++){
       
       /*Ly: Ex^p(i,j), Ex^p(i,j+1), Ex^p(i-1,j+1), Ex^p(i-1,j) */
       list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EX][i][j], -C*C*gm*bty1[j]*btx0[i]/Dx/Dy));
       list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EX][i][j+1], C*C*gm*bty1[j]*btx0[i]/Dx/Dy));
       list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EX][i-1][j+1], -C*C*gm*bty1[j]*btx0[i]/Dx/Dy));
       list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EX][i-1][j], C*C*gm*bty1[j]*btx0[i]/Dx/Dy));

        /*Ly: Ey^p(i,j), Ey^p(i+1,j), Ey^p(i-1,j) */
       list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EY][i][j],1+C*C*gm*btx0[i]/Dx/Dx*(btx1[i]+btx1[i-1]) ));
       list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EY][i+1][j], -C*C*gm*btx0[i]*btx1[i]/Dx/Dx) );
       list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EY][i-1][j],  -C*C*gm*btx0[i]*btx1[i-1]/Dx/Dx) );
       
     }
   }


   for(int i = ix; i <= ie; i++){
     for(int j = 1; j <= jh; j++){
       
       /*Lz: Ez^p(i,j), Ez^p(i,j+1), Ez^p(i,j-1), Ez^p(i+1,j), Ez^p(i-1,j) */
       list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i][j], 1 + C*C*gm*bty0[j]/Dy/Dy*(bty1[j]+bty1[j-1]) + C*C*gm*btx0[i]/Dx/Dx*(btx1[i]+btx1[i-1]) ) );
       list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i][j+1],  -C*C*gm*bty0[j]*bty1[j]/Dy/Dy));
       list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i][j-1],  -C*C*gm*bty0[j]*bty1[j-1]/Dy/Dy));
       list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i+1][j],  -C*C*gm*btx1[i]*btx0[i]/Dx/Dx));
       list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i-1][j],  -C*C*gm*btx1[i-1]*btx0[i]/Dx/Dx));
    
     }
   }





   
   /* (L+1,1)~(Nx-L,L) */
   
   for(int i = L; i < ix; i++){
     for(int j = 1; j <= L; j++){
       
       /*Lx: Ex^p(i,j), Ex^p(i,j+1), Ex^p(i,j-1) */
       list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EX][i][j] , 1+C*C*gm*bty0[j]/Dy/Dy*(bty1[j]+bty1[j-1])));
       list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EX][i][j+1], -C*C*gm*bty0[j]*bty1[j]/Dy/Dy));
       list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EX][i][j-1], -C*C*gm*bty0[j]*bty1[j-1]/Dy/Dy));

       /*Lx: Ey^p(i,j), Ey^p(i,j-1), Ey^p(i+1,j), Ey^p(i+1,j-1) */
       list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EY][i][j], -C*C*gm*bty0[j]*btx1[i]/Dx/Dy) );
       list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EY][i][j-1],  C*C*gm*bty0[j]*btx1[i]/Dx/Dy) );
       list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EY][i+1][j],  C*C*gm*bty0[j]*btx1[i]/Dx/Dy) );
       list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EY][i+1][j-1],  -C*C*gm*bty0[j]*btx1[i]/Dx/Dy) );
     }
   }

   
   for(int i = L+1; i < ix; i++){
     for(int j = 0; j < L; j++){
       /*Ly: Ex^p(i,j), Ex^p(i,j+1), Ex^p(i-1,j+1), Ex^p(i-1,j) */
       list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EX][i][j], -C*C*gm*bty1[j]*btx0[i]/Dx/Dy));
       list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EX][i][j+1], C*C*gm*bty1[j]*btx0[i]/Dx/Dy));
       list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EX][i-1][j+1], -C*C*gm*bty1[j]*btx0[i]/Dx/Dy));
       list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EX][i-1][j], C*C*gm*bty1[j]*btx0[i]/Dx/Dy));

       /*Ly: Ey^p(i,j), Ey^p(i+1,j), Ey^p(i-1,j) */
       list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EY][i][j],1+C*C*gm*btx0[i]/Dx/Dx*(btx1[i]+btx1[i-1]) ));
       list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EY][i+1][j], -C*C*gm*btx0[i]*btx1[i]/Dx/Dx) );
       list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EY][i-1][j],  -C*C*gm*btx0[i]*btx1[i-1]/Dx/Dx) );
       
     }
   }

   
   for(int i = L+1; i < ix; i++){
     for(int j = 1; j <= L; j++){
       
       /*Lz: Ez^p(i,j), Ez^p(i,j+1), Ez^p(i,j-1), Ez^p(i+1,j), Ez^p(i-1,j) */
       list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i][j], 1 + C*C*gm*bty0[j]/Dy/Dy*(bty1[j]+bty1[j-1]) + C*C*gm*btx0[i]/Dx/Dx*(btx1[i]+btx1[i-1]) ) );
       list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i][j+1],  -C*C*gm*bty0[j]*bty1[j]/Dy/Dy));
       list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i][j-1],  -C*C*gm*bty0[j]*bty1[j-1]/Dy/Dy));
       list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i+1][j],  -C*C*gm*btx1[i]*btx0[i]/Dx/Dx));
       list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i-1][j],  -C*C*gm*btx1[i-1]*btx0[i]/Dx/Dx));
      
  
     }
   }


  
 

  
   /*Ionoの範囲*/
   
   int jb {(int)(Ry/Dy)-L};
   
   /* (L+1,jh+1)~(Nx-L-1,Ny-L-1) */
   
   for(int i = L; i < ix; i++){
     for(int j = jh+1; j < jb; j++){
       
       /*Lx: Ex^p(i,j), Ex^p(i,j-1), Ex^p(i,j+1) */
       list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EX][i][j] , 1+P[0][0]+8*C*C/s/s/Dy/Dy));
       list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EX][i][j+1], -4*C*C/s/s/Dy/Dy));
       list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EX][i][j-1], -4*C*C/s/s/Dy/Dy));

         
       /*Lx: Ey^p(i,j), Ey^p(i,j-1), Ey^p(i+1,j), Ey^p(i+1,j-1) */
       list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EY][i][j], P[0][1]/4 - 4*C*C/s/s/Dx/Dy) );
       list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EY][i][j-1], P[0][1]/4 + 4*C*C/s/s/Dx/Dy) );
       list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EY][i+1][j],  P[0][1]/4 + 4*C*C/s/s/Dx/Dy) );
       list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EY][i+1][j-1],  P[0][1]/4 - 4*C*C/s/s/Dx/Dy) );

       /*Lx: Ez^p(i,j), Ez^p(i+1,j) */
       list.push_back( Eigen::Triplet <double> ( nd[EX][i][j], nd[EZ][i][j], P[0][2]/2) );
       list.push_back( Eigen::Triplet <double> ( nd[EX][i][j], nd[EZ][i+1][j], P[0][2]/2) );

     }
   }

   
   for(int i = L+1; i < ix; i++){
     for(int j = jh; j < jb; j++){
       
       /*Ly: Ex^p(i,j), Ex^p(i,j+1), Ex^p(i-1,j+1), Ex^p(i-1,j) */
       list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EX][i][j], P[1][0]/4 - 4*C*C/s/s/Dx/Dy));
       list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EX][i][j+1], P[1][0]/4 + 4*C*C/s/s/Dx/Dy));
       list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EX][i-1][j+1], P[1][0]/4 - 4*C*C/s/s/Dx/Dy));
       list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EX][i-1][j], P[1][0]/4 + 4*C*C/s/s/Dx/Dy));
       
       /*Ly: Ey^p(i,j), Ey^p(i+1,j), Ey^p(i-1,j) */
       list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EY][i][j], 1+P[1][1]+8*C*C/s/s/Dx/Dx) );
       list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EY][i+1][j], -4*C*C/s/s/Dx/Dx) );
       list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EY][i-1][j], -4*C*C/s/s/Dx/Dx) );

        /*Ly: Ez^p(i,j), Ez^p(i,j+1) */
       list.push_back( Eigen::Triplet <double> ( nd[EY][i][j], nd[EZ][i][j], P[1][2]/2) );
       list.push_back( Eigen::Triplet <double> ( nd[EY][i][j], nd[EZ][i][j+1], P[1][2]/2) );
     }
   }

   
   for(int i = L+1; i < ix; i++){
     for(int j = jh+1; j < jb; j++){
       
       /*Lz: Ex^p(i,j), Ex^p(i-1, j) */
       list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j],  nd[EX][i][j] , P[2][0]/2));
       list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j],  nd[EX][i-1][j] , P[2][0]/2));

        /*Lz: Ey^p(i,j), Ey^p(i,j-1) */
       list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j],  nd[EY][i][j] , P[2][1]/2));
       list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j],  nd[EY][i][j-1] , P[2][1]/2));

         /*Lz: Ez^p(i,j), Ez^p(i,j+1), Ez^p(i,j-1), Ez^p(i+1,j), Ez^p(i-1,j) */
       list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i][j], 1+P[2][2]+8*C*C/s/s/Dx/Dx+8*C*C/s/s/Dy/Dy) );
       list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i][j+1],  -4*C*C/s/s/Dy/Dy) );
       list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i][j-1],  -4*C*C/s/s/Dy/Dy) );
       list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i+1][j],  -4*C*C/s/s/Dx/Dx) );
       list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i-1][j],  -4*C*C/s/s/Dx/Dx) );
     }
   }
   
 



    
  /*Iono_PMLの範囲*/

    int jp {(int)(Ry/Dy)-1};
    int jq {(int)(Ry/Dy)-L};
    
    /* (1,jh+1) ~ (L,jp)*/
    
    for(int i = 0; i < L; i++){
      for(int j = jh + 1; j <= jp; j++){
	
	/*Lx: Ex^p(i,j), Ex^p(i,j+1), Ex^p(i,j-1) */
	list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EX][i][j] , 1+P[0][0]+C*C*bty0[j]/Dy/Dy*(bty1[j] + bty1[j-1])));
	list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EX][i][j+1], -C*C*bty0[j]*bty1[j]/Dy/Dy));
	list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EX][i][j-1], -C*C*bty0[j]*bty1[j-1]/Dy/Dy));

	/*Lx: Ey^p(i,j), Ey^p(i,j-1), Ey^p(i+1,j), Ey^p(i+1,j-1) */
	list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EY][i][j], P[0][1]/4 - C*C/bty0[j]*btx1[i]/Dx/Dy) );
	list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EY][i][j-1], P[0][1]/4 +  C*C/bty0[j]*btx1[i]/Dx/Dy) );
	list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EY][i+1][j],  P[0][1]/4 + C*C/bty0[j]*btx1[i]/Dx/Dy) );
	list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EY][i+1][j-1],  P[0][1]/4 - C*C/bty0[j]*btx1[i]/Dx/Dy) );

	/*Lx: Ez^p(i,j), Ez^p(i+1,j) */
	list.push_back( Eigen::Triplet <double> ( nd[EX][i][j], nd[EZ][i][j], P[0][2]/2) );
	list.push_back( Eigen::Triplet <double> ( nd[EX][i][j], nd[EZ][i+1][j], P[0][2]/2) );
      }
    }

    
    for(int i = 1; i <= L; i++){
      for(int j = jh; j <= jp; j++){
	
	/*Ly: Ex^p(i,j), Ex^p(i,j+1), Ex^p(i-1,j+1), Ex^p(i-1,j) */
	list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EX][i][j], P[1][0]/4 - C*C*btx0[i]*bty1[j]/Dx/Dy));
	list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EX][i][j+1], P[1][0]/4 + C*C*btx0[i]*bty1[j]/Dx/Dy));
	list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EX][i-1][j+1], P[1][0]/4 - C*C*btx0[i]*bty1[j]/Dx/Dy));
	list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EX][i-1][j], P[1][0]/4 + C*C*btx0[i]*bty1[j]/Dx/Dy));

	/*Ly: Ey^p(i,j), Ey^p(i+1,j), Ey^p(i-1,j) */
	list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EY][i][j], 1+P[1][1]+C*C/btx0[i]/Dx/Dx*(btx1[i]+btx1[i-1])) );
	list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EY][i+1][j], -C*C/btx0[i]/Dx/Dx*btx1[i] ));
	list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EY][i-1][j], -C*C/btx0[i]/Dx/Dx*btx1[i-1] ));

	/*Ly: Ez^p(i,j), Ez^p(i,j+1) */
	list.push_back( Eigen::Triplet <double> ( nd[EY][i][j], nd[EZ][i][j], P[1][2]/2) );
	list.push_back( Eigen::Triplet <double> ( nd[EY][i][j], nd[EZ][i][j+1], P[1][2]/2) );
      }
    }

    
    for(int i = 1; i <= L; i++){
      for(int j = jh + 1; j <= jp; j++){
	/*Lz: Ex^p(i,j) */
        list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j],  nd[EX][i][j] , P[2][0]/2));
        list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j],  nd[EX][i-1][j] , P[2][0]/2));

	/*Lz: Ey^p(i,j), Ey^p(i,j-1) */
	list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j],  nd[EY][i][j] , P[2][1]/2));
	list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j],  nd[EY][i][j-1] , P[2][1]/2));

	/*Lz: Ez^p(i,j), Ez^p(i,j+1), Ez^p(i,j-1), Ez^p(i+1,j), Ez^p(i-1,j) */
	list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i][j], 1+P[2][2]+C*C*btx0[i]/Dx/Dx*(btx1[i] + btx1[i-1])+C*C*bty0[j]/Dy/Dy*(bty1[j] + bty1[j-1])) );
	list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i][j+1],  -C*C*bty0[j]*bty1[j]/Dy/Dy) );
	list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i][j-1],  -C*C*bty0[j]*bty1[j-1]/Dy/Dy) );
	list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i+1][j],   -C*C*btx0[i]*btx1[i]/Dx/Dx) );
	list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i-1][j],   -C*C*btx0[i]*btx1[i-1]/Dx/Dx) );
      }
    }
    


    

    
    /* (Nx-L,jh+1) ~ (ie,jp)*/
    
    for(int i = ix; i <= ie; i++){
      for(int j = jh + 1; j <= jp; j++){
	
	/*Lx: Ex^p(i,j), Ex^p(i,j+1), Ex^p(i,j-1) */
	list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EX][i][j] , 1+P[0][0]+C*C*bty0[j]/Dy/Dy*(bty1[j] + bty1[j-1])));
	list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EX][i][j+1], -C*C*bty0[j]*bty1[j]/Dy/Dy));
	list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EX][i][j-1], -C*C*bty0[j]*bty1[j-1]/Dy/Dy));

	/*Lx: Ey^p(i,j), Ey^p(i,j-1), Ey^p(i+1,j), Ey^p(i+1,j-1) */
	list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EY][i][j], P[0][1]/4 - C*C/bty0[j]*btx1[i]/Dx/Dy) );
	list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EY][i][j-1], P[0][1]/4 +  C*C/bty0[j]*btx1[i]/Dx/Dy) );
	list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EY][i+1][j],  P[0][1]/4 + C*C/bty0[j]*btx1[i]/Dx/Dy) );
	list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EY][i+1][j-1],  P[0][1]/4 - C*C/bty0[j]*btx1[i]/Dx/Dy) );

	/*Lx: Ez^p(i,j), Ez^p(i+1,j) */
	list.push_back( Eigen::Triplet <double> ( nd[EX][i][j], nd[EZ][i][j], P[0][2]/2) );
	list.push_back( Eigen::Triplet <double> ( nd[EX][i][j], nd[EZ][i+1][j], P[0][2]/2) );
      }
    }

    
    for(int i = ix; i <= ie; i++){
      for(int j = jh; j <= jp; j++){	 
	/*Ly: Ex^p(i,j), Ex^p(i,j+1), Ex^p(i-1,j+1), Ex^p(i-1,j) */
	list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EX][i][j], P[1][0]/4 - C*C*btx0[i]*bty1[j]/Dx/Dy));
	list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EX][i][j+1], P[1][0]/4 + C*C*btx0[i]*bty1[j]/Dx/Dy));
	list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EX][i-1][j+1], P[1][0]/4 - C*C*btx0[i]*bty1[j]/Dx/Dy));
	list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EX][i-1][j], P[1][0]/4 + C*C*btx0[i]*bty1[j]/Dx/Dy));

	/*Ly: Ey^p(i,j), Ey^p(i+1,j), Ey^p(i-1,j) */
	list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EY][i][j], 1+P[1][1]+C*C/btx0[i]/Dx/Dx*(btx1[i]+btx1[i-1])) );
	list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EY][i+1][j], -C*C/btx0[i]/Dx/Dx*btx1[i] ));
	list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EY][i-1][j], -C*C/btx0[i]/Dx/Dx*btx1[i-1] ));

	/*Ly: Ez^p(i,j), Ez^p(i,j+1) */
	list.push_back( Eigen::Triplet <double> ( nd[EY][i][j], nd[EZ][i][j], P[1][2]/2) );
	list.push_back( Eigen::Triplet <double> ( nd[EY][i][j], nd[EZ][i][j+1], P[1][2]/2) );
      }
    }

    
    for(int i = ix; i <= ie; i++){
      for(int j = jh + 1; j <= jp; j++){
	/*Lz: Ex^p(i,j) */
        list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j],  nd[EX][i][j] , P[2][0]/2));
        list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j],  nd[EX][i-1][j] , P[2][0]/2));

	/*Lz: Ey^p(i,j), Ey^p(i,j-1) */
	list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j],  nd[EY][i][j] , P[2][1]/2));
	list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j],  nd[EY][i][j-1] , P[2][1]/2));

	/*Lz: Ez^p(i,j), Ez^p(i,j+1), Ez^p(i,j-1), Ez^p(i+1,j), Ez^p(i-1,j) */
	list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i][j], 1+P[2][2]+C*C*btx0[i]/Dx/Dx*(btx1[i] + btx1[i-1])+C*C*bty0[j]/Dy/Dy*(bty1[j] + bty1[j-1])) );
	list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i][j+1],  -C*C*bty0[j]*bty1[j]/Dy/Dy) );
	list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i][j-1],  -C*C*bty0[j]*bty1[j-1]/Dy/Dy) );
	list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i+1][j],   -C*C*btx0[i]*btx1[i]/Dx/Dx) );
	list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i-1][j],   -C*C*btx0[i]*btx1[i-1]/Dx/Dx) );
      }
    }
    



    
    
    /* (L+1,jq) ~ (ix,jp)*/
    
    for(int i = L; i < ix; i++){
      for(int j = jq; j <= jp; j++){
	
	/*Lx: Ex^p(i,j), Ex^p(i,j+1), Ex^p(i,j-1) */
	list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EX][i][j] , 1+P[0][0]+C*C*bty0[j]/Dy/Dy*(bty1[j] + bty1[j-1])));
	list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EX][i][j+1], -C*C*bty0[j]*bty1[j]/Dy/Dy));
	list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EX][i][j-1], -C*C*bty0[j]*bty1[j-1]/Dy/Dy));

	/*Lx: Ey^p(i,j), Ey^p(i,j-1), Ey^p(i+1,j), Ey^p(i+1,j-1) */
	list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EY][i][j], P[0][1]/4 - C*C/bty0[j]*btx1[i]/Dx/Dy) );
	list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EY][i][j-1], P[0][1]/4 +  C*C/bty0[j]*btx1[i]/Dx/Dy) );
	list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EY][i+1][j],  P[0][1]/4 + C*C/bty0[j]*btx1[i]/Dx/Dy) );
	list.push_back( Eigen::Triplet <double> ( nd[EX][i][j],  nd[EY][i+1][j-1],  P[0][1]/4 - C*C/bty0[j]*btx1[i]/Dx/Dy) );

	/*Lx: Ez^p(i,j), Ez^p(i+1,j) */
	list.push_back( Eigen::Triplet <double> ( nd[EX][i][j], nd[EZ][i][j], P[0][2]/2) );
	list.push_back( Eigen::Triplet <double> ( nd[EX][i][j], nd[EZ][i+1][j], P[0][2]/2) );
      }
    }

    
    for(int i = L+1; i < ix; i++){
      for(int j = jq; j <= jp; j++){
	/*Ly: Ex^p(i,j), Ex^p(i,j+1), Ex^p(i-1,j+1), Ex^p(i-1,j) */
	list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EX][i][j], P[1][0]/4 - C*C*btx0[i]*bty1[j]/Dx/Dy));
	list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EX][i][j+1], P[1][0]/4 + C*C*btx0[i]*bty1[j]/Dx/Dy));
	list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EX][i-1][j+1], P[1][0]/4 - C*C*btx0[i]*bty1[j]/Dx/Dy));
	list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EX][i-1][j], P[1][0]/4 + C*C*btx0[i]*bty1[j]/Dx/Dy));
	
	/*Ly: Ey^p(i,j), Ey^p(i+1,j), Ey^p(i-1,j) */
	list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EY][i][j], 1+P[1][1]+C*C/btx0[i]/Dx/Dx*(btx1[i]+btx1[i-1])) );
	list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EY][i+1][j], -C*C/btx0[i]/Dx/Dx*btx1[i] ));
	list.push_back( Eigen::Triplet <double> ( nd[EY][i][j],  nd[EY][i-1][j], -C*C/btx0[i]/Dx/Dx*btx1[i-1] ));

	/*Ly: Ez^p(i,j), Ez^p(i,j+1) */
	list.push_back( Eigen::Triplet <double> ( nd[EY][i][j], nd[EZ][i][j], P[1][2]/2) );
	list.push_back( Eigen::Triplet <double> ( nd[EY][i][j], nd[EZ][i][j+1], P[1][2]/2) );
      }
    }
    
    for(int i = L+1; i <= ix; i++){
      for(int j = jq; j <= jp; j++){
	/*Lz: Ex^p(i,j) */
        list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j],  nd[EX][i][j] , P[2][0]/2));
        list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j],  nd[EX][i-1][j] , P[2][0]/2));

	/*Lz: Ey^p(i,j), Ey^p(i,j-1) */
	list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j],  nd[EY][i][j] , P[2][1]/2));
	list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j],  nd[EY][i][j-1] , P[2][1]/2));

	/*Lz: Ez^p(i,j), Ez^p(i,j+1), Ez^p(i,j-1), Ez^p(i+1,j), Ez^p(i-1,j) */
	list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i][j], 1+P[2][2]+C*C*btx0[i]/Dx/Dx*(btx1[i] + btx1[i-1])+C*C*bty0[j]/Dy/Dy*(bty1[j] + bty1[j-1])) );
	list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i][j+1],  -C*C*bty0[j]*bty1[j]/Dy/Dy) );
	list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i][j-1],  -C*C*bty0[j]*bty1[j-1]/Dy/Dy) );
	list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i+1][j],   -C*C*btx0[i]*btx1[i]/Dx/Dx) );
	list.push_back( Eigen::Triplet <double> ( nd[EZ][i][j], nd[EZ][i-1][j],   -C*C*btx0[i]*btx1[i-1]/Dx/Dx) );
      }
    }
    

    /*
    list.push_back( Eigen::Triplet <double> (nd[EX][0][0], nd[EX][0][0], 1.0) );
    list.push_back( Eigen::Triplet <double> (nd[EX][Nx-1][Ny], nd[EX][Nx-1][Ny], 1.0) );
    list.push_back( Eigen::Triplet <double> (nd[EY][0][0], nd[EY][0][0], 1.0) );
    list.push_back( Eigen::Triplet <double> (nd[EY][Nx][Ny-1], nd[EY][Nx][Ny-1], 1.0) );
    list.push_back( Eigen::Triplet <double> (nd[EZ][0][0], nd[EZ][0][0], 1.0) );
    list.push_back( Eigen::Triplet <double> (nd[EZ][Nx][Ny], nd[EZ][Nx][Ny], 1.0) );
    */

    
    
    for(int i = 0; i < Nx; i++){
      list.push_back( Eigen::Triplet <double> (nd[EX][i][0], nd[EX][i][0], 1.0) );
      list.push_back( Eigen::Triplet <double> (nd[EX][i][Ny], nd[EX][i][Ny], 1.0) );
    }
    
    
    for(int j = 0; j < Ny; j++){
      list.push_back( Eigen::Triplet <double> (nd[EY][0][j], nd[EY][0][j], 1.0) );
      list.push_back( Eigen::Triplet <double> (nd[EY][Nx][j], nd[EY][Nx][j], 1.0) );
    }

    for(int i = 0; i <= Nx; i++){
      list.push_back( Eigen::Triplet <double> (nd[EZ][i][0], nd[EZ][i][0], 1.0) );
      list.push_back( Eigen::Triplet <double> (nd[EZ][i][Ny], nd[EZ][i][Ny], 1.0) );
    }
    
    for(int j = 0; j <= Ny; j++){
      list.push_back( Eigen::Triplet <double> (nd[EZ][0][j], nd[EZ][0][j], 1.0) );
      list.push_back( Eigen::Triplet <double> (nd[EZ][Nx][j], nd[EZ][Nx][j], 1.0) );
    }
    
   

    
    
    Eigen::SparseMatrix <double> A( Nx*(Ny+1) + (Nx+1)*Ny + (Nx+1)*(Ny+1),  Nx*(Ny+1) + (Nx+1)*Ny + (Nx+1)*(Ny+1));
    
    A.setFromTriplets (list.begin(), list.end() );
    A.makeCompressed();
    
    S.analyzePattern(A);
    S.factorize(A);


    delete1d(btx0);
    delete1d(btx1);
    delete1d(bty0);
    delete1d(bty1);

}

