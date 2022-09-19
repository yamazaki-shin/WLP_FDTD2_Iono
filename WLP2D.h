#include <eigen3/Eigen/Sparse>
#include <cmath>

constexpr double C0 { 3.0e8 }; /* 真空中の光速 */
constexpr double MU0 { 4.0 * M_PI * 1.0e-7 }; /* 真空の透磁率 */
constexpr double EPS0 { 1.0 / MU0 / C0 / C0 }; /* 真空の誘電率 */
constexpr double SIG0 { 0.0 }; /* 導電率 */

constexpr double Ne { 1.0e8 }; /* 電子密度 */
constexpr double B0 { 5.0e-5 };
constexpr double m { 9.109e-31 }; /* 電子の質量 */
constexpr double e { 1.602e-19 }; /* 電子の電荷 */
constexpr double wp {sqrt(e*e*Ne/EPS0/m)}; /*プラズマ周波数*/
constexpr double wc { e*B0/m };  /*サイクロトロン周波数*/ 

constexpr int Nx { 100 };
constexpr int Ny { 100 };
constexpr double Rx { 100.0e3}; /* 計算する範囲 */
constexpr double Ry { 100.0e3}; /* 計算する範囲 */
constexpr double Dx { Rx / double(Nx) };
constexpr double Dy { Ry / double(Ny) };


constexpr double Dt { 0.99/(C0*sqrt((1.0/Dx)*(1.0/Dx) + (1.0/Dy)*(1.0/Dy)))};
constexpr double St { 10.0* Dt}; /*σt*/
constexpr double t0 { 6.0* St};


constexpr int P {50}; /*次数*/
constexpr int M {1000}; 
constexpr int N {200};
constexpr int L {8}; /*PMLの範囲*/
constexpr double dt  {2*t0/M};
constexpr double Ez_dt  {4*t0/N};
constexpr double s {1.0e4}; /*スケールパラメ➖タ*/
constexpr double z {70}; /* 高度 */


constexpr int EX {0};  /*ノード番号*/
constexpr int EY {1};
constexpr int EZ {2};


constexpr double Sx {50e3};
constexpr double Sy {20e3+L*Dy};
constexpr double hi {64e3};


constexpr double a {1/(EPS0*s/2+SIG0)};
constexpr double b {2/MU0/s};
constexpr double C {1/MU0/s};
constexpr double gm {s*EPS0/(1+2*SIG0)};


double ***memory_allocate3d(int M, int N, int L, double initial_value);
int ***memory_allocate3d_int(int M, int N, int L,  int initial_value);
double **memory_allocate2d(int M, int N, double initial_value);
double *memory_allocate1d(int M, double initial_value);
void delete3d(double ***v, int M, int N);
void delete3d_int(int ***v, int M, int N);
void delete2d(double **v, int M);
void delete1d(double *v);

void cal_Hxp(double **Hxp, double **SUM_Hxp,double **Ezp);
void cal_Hyp(double **Hyp, double **SUM_Hyp,double **Ezp);
void cal_Hzp(double **Hzp,double **Hzxp, double **Hzyp, double **SUM_Hzp,double **Exp, double **Eyp);
void cal_Hzxp(double **Hzxp, double **SUM_Hzxp, double **Eyp);
void cal_Hzyp(double **Hzyp, double **SUM_Hzyp, double **Exp);

void cal_Fxp(double **Fxp, double **SUM_Fxp,double **Hzp);
void cal_Fyp(double **Fyp, double **SUM_Fyp,double **Hzp);
void cal_Fzxp(double **Fzxp, double **SUM_Fzxp, double **Hyp);
void cal_Fzyp(double **Fzyp, double **SUM_Fzyp, double **Hxp);

void cal_Jmxp(double **Jmxp, double **SUM_Jmxp, double **SUM_Jmyp, double **SUM_Jmzp, double **Exp, double **Eyp, double **Ezp);
void cal_Jmyp(double **Jmyp, double **SUM_Jmxp, double **SUM_Jmyp, double **SUM_Jmzp, double **Exp, double **Eyp, double **Ezp);
void cal_Jmzp(double **Jmzp, double **SUM_Jmxp, double **SUM_Jmyp, double **SUM_Jmzp, double **Exp, double **Eyp, double **Ezp);


void cal_RHS_vector(Eigen::VectorXd &b, double **SUM_Hxp, double **SUM_Hyp, double **SUM_Hzp, double **SUM_Hzxp, double **SUM_Hzyp, double **SUM_Exp, double **SUM_Eyp, double **SUM_Ezp,double **SUM_Fxp, double **SUM_Fyp, double **SUM_Fzxp, double **SUM_Fzyp, double **SUM_Jmxp, double **SUM_Jmyp, double **SUM_Jmzp, double *Jzp, int q, int ***nd);
void sum(double **SUM, double **f, int n, int m);
double Jz(double t);
void cal_Jz_coefficient(double *Jzp);
void Fai(double t, double *F);
void cal_Ez(double ***Ez,double **Ezp,int n);
void compose_coef_matrix(Eigen::SparseLU < Eigen::SparseMatrix <double>, Eigen::COLAMDOrdering <int> > &S,  int ***nd, double **P);
void cal_Exp(double **Ezp, int ***nd, Eigen::VectorXd &b, Eigen::SparseLU < Eigen::SparseMatrix <double>, Eigen::COLAMDOrdering <int> > &S);
void cal_Eyp(double **Ezp, int ***nd, Eigen::VectorXd &b, Eigen::SparseLU < Eigen::SparseMatrix <double>, Eigen::COLAMDOrdering <int> > &S);
void cal_Ezp(double **Ezp, int ***nd, Eigen::VectorXd &b, Eigen::SparseLU < Eigen::SparseMatrix <double>, Eigen::COLAMDOrdering <int> > &S);
void cal_nd(int ***nd);
void PML_cond(double *btx0, double *btx1, double *bty0, double *bty1);

void cal_Nyu(double v, int z);
void cal_P(double **P);
