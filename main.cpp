#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <eigen3/Eigen/Sparse>
#include "WLP2D.h"


int main(void){
  
  double ***Ex = memory_allocate3d(Nx+1, Ny+1, N+1, 0.0); /*Ex(x, y, t)*の3次元配列*/
  double ***Ey = memory_allocate3d(Nx+1, Ny+1, N+1, 0.0); /*Ey(x, y, t)*の3次元配列*/
  double ***Ez = memory_allocate3d(Nx+1, Ny+1, N+1, 0.0); /*Ez(x, y, t)*の3次元配列*/
  
  int ***nd = memory_allocate3d_int(3, Nx+1, Ny+1, -1); /*ノード番号の3次元配列*/
  
  double **Exp = memory_allocate2d(Nx, Ny+1, 0.0); /*Ex^p(x, y)*の2次元配列*/
  double **Eyp = memory_allocate2d(Nx+1, Ny, 0.0); /*Ey^p(x, y)*の2次元配列*/
  double **Ezp = memory_allocate2d(Nx+1, Ny+1, 0.0);  /*Ez^p(x, y)*の2次元配列*/
  double **Exp_S = memory_allocate2d(Nx+1, Ny+1, 0.0); /*Ex^p-1_SUMの2次元配列*/
  double **Eyp_S = memory_allocate2d(Nx+1, Ny+1, 0.0); /*Ey^p-1_SUMの2次元配列*/
  double **Ezp_S = memory_allocate2d(Nx+1, Ny+1, 0.0); /*Ez^p-1_SUMの2次元配列*/
  
  double **Hxp =  memory_allocate2d(Nx+1, Ny, 0.0);  /*Hx^p(x, y)*の2次元配列*/
  double **Hyp =  memory_allocate2d(Nx, Ny+1, 0.0);  /*Hy^p(x, y)*の2次元配列*/
  double **Hzp =  memory_allocate2d(Nx, Ny, 0.0);  /*Hz^p(x, y)*の2次元配列*/
  double **Hzxp =  memory_allocate2d(Nx, Ny, 0.0);  /*Hzx^p(x, y)*の2次元配列*/
  double **Hzyp =  memory_allocate2d(Nx, Ny, 0.0);  /*Hzy^p(x, y)*の2次元配列*/
  double **Hxp_S =  memory_allocate2d(Nx+1, Ny, 0.0);  /*Hx^p-1_SUM(x, y)*の2次元配列*/
  double **Hyp_S =  memory_allocate2d(Nx, Ny+1, 0.0);  /*Hy^p-1_SUM(x, y)*の2次元配列*/
  double **Hzp_S =  memory_allocate2d(Nx, Ny, 0.0);  /*Hz^p-1_SUM(x, y)*の2次元配列*/
  double **Hzxp_S =  memory_allocate2d(Nx, Ny, 0.0);  /*Hzx^p-1_SUM(x, y)*の2次元配列*/
  double **Hzyp_S =  memory_allocate2d(Nx, Ny, 0.0);  /*Hzy^p-1_SUM(x, y)*の2次元配列*/
  
  double **Fxp = memory_allocate2d(Nx, Ny+1, 0.0); /*Fx^p(x, y)*の2次元配列*/
  double **Fyp = memory_allocate2d(Nx+1, Ny, 0.0); /*Fy^p(x, y)*の2次元配列*/
  double **Fzxp = memory_allocate2d(Nx+1, Ny+1, 0.0);  /*Fzx^p(x, y)*の2次元配列*/
  double **Fzyp = memory_allocate2d(Nx+1, Ny+1, 0.0); /*Fzy^p(x, y)*の2次元配列*/
  double **Fxp_S =  memory_allocate2d(Nx+1, Ny+1, 0.0);  /*Fx^p-1_SUM(x, y)*の2次元配列*/
  double **Fyp_S =  memory_allocate2d(Nx+1, Ny+1, 0.0);  /*Fy^p-1_SUM(x, y)*の2次元配列*/
  double **Fzxp_S =  memory_allocate2d(Nx+1, Ny+1, 0.0);  /*Fzx^p-1_SUM(x, y)*の2次元配列*/
  double **Fzyp_S =  memory_allocate2d(Nx+1, Ny+1, 0.0);  /*Fzy^p-1_SUM(x, y)*の2次元配列*/
  
  double **Pm =  memory_allocate2d(3, 3, 0.0);  /*P(x, y)*の2次元配列*/
  
  double **Jmxp =  memory_allocate2d(Nx+1, Ny+1, 0.0);  /*Jmx^p(x, y)*の2次元配列*/
  double **Jmyp =  memory_allocate2d(Nx+1, Ny+1, 0.0);  /*Jmy^p(x, y)*の2次元配列*/
  double **Jmzp =  memory_allocate2d(Nx+1, Ny+1, 0.0);  /*Jmz^p(x, y)*の2次元配列*/
  double **Jmxp_S =  memory_allocate2d(Nx+1, Ny+1, 0.0);  /*Jmx^p-1_SUM(x, y)*の2次元配列*/
  double **Jmyp_S =  memory_allocate2d(Nx+1, Ny+1, 0.0);  /*Jmy^p-1_SUM(x, y)*の2次元配列*/
  double **Jmzp_S =  memory_allocate2d(Nx+1, Ny+1, 0.0);  /*Jmz^p-1_SUM(x, y)*の2次元配列*/
  
  double *Jzp =  memory_allocate1d(P+1, 0.0); /*各時刻の電流密度を確保する配列*/

  cal_nd(nd);
  cal_P(Pm);  /*プラズマ電流*/
  // for(int i = 0; i <= 2; i++){
  //   for(int j = 0; j <= Ny; j++){
  // for(int n = 0; n <=2; n++){
  //	std::cout << n <<" " << i << " " << j << " " << nd[n][i][j] << std::endl;
  //  } 
  //  }
  //  }

   cal_Jz_coefficient(Jzp); /*展開係数Jzpを計算*/

   Eigen::VectorXd b = Eigen::VectorXd::Zero( Nx*(Ny+1) + (Nx+1)*Ny + (Nx+1)*(Ny+1));
   Eigen::SparseLU < Eigen::SparseMatrix <double>, Eigen::COLAMDOrdering <int> >  S;
   
   compose_coef_matrix(S,nd,Pm);    /*係数行列を生成*/
   
   for(int q = 0; q <= P; q++){
     
     cal_RHS_vector(b, Hxp_S, Hyp_S, Hzp_S, Hzxp_S, Hzyp_S, Exp_S, Eyp_S, Ezp_S, Fxp_S, Fyp_S, Fzxp_S, Fzyp_S, Jmxp_S, Jmyp_S, Jmzp_S, Jzp, q, nd); /*右辺ベクトル計算*/
     // std::cout << "a" << std::endl;
     /*Exp, Eyp, Ezpを解く*/
     //   cal_Exp(Exp, nd, b, S);
     //  cal_Eyp(Eyp, nd, b, S);
     //cal_Ezp(Ezp, nd, b, S);
     // std::cout << "b" << std::endl;
     
     
     /*Ezを計算*/
     cal_Ez(Ez,Ezp, q);
     //std::cout << "Ez" << Ez[50][50][50] << std::endl;    
     /*Hxp,Hyp,Hzp,Hzxp,Hzyp,Fxp,Fyp,Fzxp,Fzypを計算*/
     cal_Hxp(Hxp, Hxp_S, Ezp);
     cal_Hyp(Hyp, Hyp_S, Ezp);
     cal_Hzp(Hzp, Hzxp, Hzyp, Hzp_S, Exp, Eyp);
     cal_Hzxp(Hzxp, Hzxp_S, Eyp);
     cal_Hzyp(Hzyp, Hzyp_S, Exp);
     cal_Fxp(Fxp, Fxp_S, Hzp);
     cal_Fyp(Fyp, Fyp_S, Hzp);
     cal_Fzxp(Fzxp, Fzxp_S, Hyp);
     cal_Fzyp(Fzyp, Fzyp_S, Hxp);
     cal_Jmxp(Jmxp, Jmxp_S, Jmyp_S, Jmzp_S, Exp, Eyp, Ezp);
     cal_Jmyp(Jmyp, Jmxp_S, Jmyp_S, Jmzp_S, Exp, Eyp, Ezp);
     cal_Jmzp(Jmzp, Jmxp_S, Jmyp_S, Jmzp_S, Exp, Eyp, Ezp);
     

    /*それぞれsumを計算*/
    sum(Exp_S, Exp,Nx,Ny+1); /*Ez^p-1_sumを計算*/
    sum(Eyp_S, Eyp,Nx+1,Ny); /*Ez^p-1_sumを計算*/
    sum(Ezp_S, Ezp,Nx+1,Ny+1); /*Ez^p-1_sumを計算*/
    sum(Hxp_S, Hxp,Nx+1,Ny); /*Hx^p-1_sumを計算*/
    sum(Hyp_S, Hyp, Nx,Ny+1); /*Hy^p-1_sumを計算*/
    sum(Hzp_S, Hzp, Nx,Ny); /*Hz^p-1_sumを計算*/
    sum(Hzxp_S, Hzxp,Nx,Ny); /*Hzx^p-1_sumを計算*/
    sum(Hzyp_S, Hzyp, Nx,Ny); /*Hzy^p-1_sumを計算*/
    sum(Jmxp_S, Jmxp,Nx+1,Ny+1); /*Jmx^p-1_sumを計算*/
    sum(Jmyp_S, Jmyp,Nx+1,Ny+1); /*Jmy^p-1_sumを計算*/
    sum(Jmzp_S, Jmzp,Nx+1,Ny+1); /*Jmz^p-1_sumを計算*/
    sum(Fxp_S, Fxp,Nx,Ny+1); /*Fx^p-1_sumを計算*/
    sum(Fyp_S, Fyp,Nx+1,Ny); /*Fy^p-1_sumを計算*/
    sum(Fzxp_S, Fzxp,Nx+1,Ny+1); /*Fzx^p-1_sumを計算*/
    sum(Fzyp_S, Fzyp,Nx+1,Ny+1); /*Fzy^p-1_sumを計算*/

    if(q % 10 == 0){
      std::cout << q << std::endl;
    }
    
  }
  
  
  for(int n = 0; n <= N ; n++){
    std::string filename = "Ez_WLP2D_Iono_" + std::to_string(n) + ".dat";
    std::ofstream ofs( filename.c_str() );
    for(int i = 0; i <= Nx; i++){
      for(int j = 0; j <= Ny; j++){
	double x = i*Dx;
	double y = j*Dy;
	ofs << x << " " << y << " " << Ez[i][j][n] << std::endl;
      }
      ofs << std::endl;
    }
    ofs.close();
  }
  
  
  
  delete3d(Ex, Nx+1,Ny+1);
  delete3d(Ey, Nx+1,Ny+1);
  delete3d(Ez, Nx+1,Ny+1);
  delete3d_int(nd,3,Nx+1);
  delete2d(Exp, Nx);
  delete2d(Eyp, Nx+1);
  delete2d(Ezp, Nx+1);  
  delete2d(Exp_S, Nx);
  delete2d(Eyp_S, Nx+1);
  delete2d(Ezp_S, Nx+1);
  delete2d(Hxp, Nx+1);
  delete2d(Hyp, Nx);
  delete2d(Hzp, Nx);
  delete2d(Hzxp, Nx);
  delete2d(Hzyp, Nx);
  delete2d(Hxp_S, Nx+1);
  delete2d(Hyp_S, Nx);
  delete2d(Hzp_S, Nx);
  delete2d(Hzxp_S, Nx);
  delete2d(Hzyp_S, Nx);
  delete1d(Jzp);
  
}

