#define _CRT_SECURE_NO_WARNINGS

#pragma once
#include <iostream>
#include <vector>
#include "const_FDTD.h"
#include "Output.h"
#include<string>
#include<omp.h>
#include<iomanip>

using namespace std;

struct PML
{
	int nx0, ny0, nz0, nx1, ny1, nz1;
	vector<vector<vector<double>>> Exy, Exz;
	vector<vector<vector<double>>> Eyz, Eyx;
	vector<vector<vector<double>>> Ezx, Ezy;
	vector<vector<vector<double>>> Hxy, Hxz;
	vector<vector<vector<double>>> Hyz, Hyx;
	vector<vector<vector<double>>> Hzx, Hzy;
	vector<vector<vector<double>>> aexy, aexz;
	vector<vector<vector<double>>> aeyx, aeyz;
	vector<vector<vector<double>>> aezx, aezy;
	vector<vector<vector<double>>> bexy, bexz;
	vector<vector<vector<double>>> beyx, beyz;
	vector<vector<vector<double>>> bezx, bezy;
	vector<vector<vector<double>>> amxy, amxz;
	vector<vector<vector<double>>> amyx, amyz;
	vector<vector<vector<double>>> amzx, amzy;
	vector<vector<vector<double>>> bmxy, bmxz;
	vector<vector<vector<double>>> bmyx, bmyz;
	vector<vector<vector<double>>> bmzx, bmzy;
};


class FDTD
{
private:
	//解析領域中の電磁界配列
	vector< vector< vector<double> > > Ex, Ey, Ez;
	vector< vector< vector<double> > > Hx, Hy, Hz;
	//PML媒質中の電磁界用の配列
	vector< vector< vector< vector<double> > > > Exyx, Exzx, Eyxx, Eyzx, Ezxx, Ezyx;
	vector< vector< vector< vector<double> > > > Exyy, Exzy, Eyxy, Eyzy, Ezxy, Ezyy;
	vector< vector< vector< vector<double> > > > Exyz, Exzz, Eyxz, Eyzz, Ezxz, Ezyz;
	vector< vector< vector< vector<double> > > > Hxyx, Hxzx, Hyxx, Hyzx, Hzxx, Hzyx;
	vector< vector< vector< vector<double> > > > Hxyy, Hxzy, Hyxy, Hyzy, Hzxy, Hzyy;
	vector< vector< vector< vector<double> > > > Hxyz, Hxzz, Hyxz, Hyzz, Hzxz, Hzyz;
	//媒質記録用の配列
	vector< vector< vector<int> > > IDE;

	//媒質の物性値
	vector<double> eps, sgmE, mu, sgmM;
	//電磁界計算用の係数
	vector<double> aex, aey, aez;
	vector<double> bexy, bexz;
	vector<double> beyx, beyz;
	vector<double> bezx, bezy;
	vector<double> amx, amy, amz;
	vector<double> bmxy, bmxz;
	vector<double> bmyx, bmyz;
	vector<double> bmzx, bmzy;
	//PML計算用の構造体（クラス）
	vector<PML> pml;

public:
	FDTD();
	void DeterminMatrixSize();
	void InitMatrix();
	void Modeling();
	void CalcCoefficient();
	void CalcCE(vector<double> &ae, vector<double> &be1, vector<double> &be2 ,double d1,double d2);
	void CalcCM(vector<double> &am, vector<double> &bm1, vector<double> &bm2, double d1, double d2);
	void CalcPMLCECM(PML *pml, int nx0, int nx1, int ny0, int ny1, int nz0, int nz1);
	void SourceE(double t);
	void SourceH(double t);
	void CalcEField();
	void CalcHField();
	void CalcPMLEField(PML *pml);
	void CalcPMLHField(PML *pml);
	void StartRepetition();

};

