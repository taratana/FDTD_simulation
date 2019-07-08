#include <math.h>

//物理定数
static const double PI = 3.14159265358979323846;
static const double EPSILON0 = 8.854e-12;
static const double MU0 = PI * 4.0e-7;
static const double C = 2.997956380e8;

//セルサイズ [m]
static const double DX = 12.5e3;	//λ/12
static const double DY = 12.5e3;
static const double DZ = 12.5e3;

//もろもろ
static const double DT = 0.99/(C*sqrt(1.0/(DX*DX)+ 1.0 / (DY*DY)+ 1.0 / (DZ*DZ)));
static const int NT = 300;
static const int NUMBER_OF_MATERIAL = 2;
static const int M = 4;			//PMLの導電率の係数
static const int L_PML = 16;	//PMLの層数
static const double REF0 = -120; //PML媒質からの反射の大きさ[dB]

//(PMLを含めない)領域分割数
static const int NXX = 120;	//解析領域に10波長収めたい
static const int NYY = 120;
static const int NZZ = 120;

//全領域(PMLを含めた)の分割数
static const int NX = NXX + 2 * L_PML;
static const int NY = NYY + 2 * L_PML;
static const int NZ = NZZ + 2 * L_PML;

//水の電気的物性値
static const double WATER_EPS = 80.36;
static const double WATER_SIGMA_E = 0.01;
static const double WATER_MU = 1.0;
static const double WATER_SIGMA_M = 0;

//波源の設定
static const double E_WAVE_AMPLITUDE = 1.0;
static const double H_WAVE_AMPLITUDE = 1.0 / sqrt(MU0 / EPSILON0);
static const double WAVE_FREQUENCY = 2000;	//[Hz]


