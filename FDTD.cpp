#include "FDTD.h"



FDTD::FDTD() {

}

void FDTD::DeterminMatrixSize() {
	Ex.resize(NX + 1, vector< vector<double> >(NY + 1, vector<double>(NZ + 1)));
	Ey.resize(NX + 1, vector< vector<double> >(NY + 1, vector<double>(NZ + 1)));
	Ez.resize(NX + 1, vector< vector<double> >(NY + 1, vector<double>(NZ + 1)));
	Hx.resize(NX + 1, vector< vector<double> >(NY + 1, vector<double>(NZ + 1)));
	Hy.resize(NX + 1, vector< vector<double> >(NY + 1, vector<double>(NZ + 1)));
	Hz.resize(NX + 1, vector< vector<double> >(NY + 1, vector<double>(NZ + 1)));

	IDE.resize(NX + 1, vector< vector<int> >(NY + 1, vector<int>(NZ + 1)));

	eps.resize(NUMBER_OF_MATERIAL);
	sgmE.resize(NUMBER_OF_MATERIAL);
	mu.resize(NUMBER_OF_MATERIAL);
	sgmM.resize(NUMBER_OF_MATERIAL);
	aex.resize(NUMBER_OF_MATERIAL); aey.resize(NUMBER_OF_MATERIAL); aez.resize(NUMBER_OF_MATERIAL);
	bexy.resize(NUMBER_OF_MATERIAL); bexz.resize(NUMBER_OF_MATERIAL);
	beyx.resize(NUMBER_OF_MATERIAL); beyz.resize(NUMBER_OF_MATERIAL);
	bezx.resize(NUMBER_OF_MATERIAL); bezy.resize(NUMBER_OF_MATERIAL);
	amx.resize(NUMBER_OF_MATERIAL);  amy.resize(NUMBER_OF_MATERIAL);  amz.resize(NUMBER_OF_MATERIAL);
	bmxy.resize(NUMBER_OF_MATERIAL); bmxz.resize(NUMBER_OF_MATERIAL);
	bmyx.resize(NUMBER_OF_MATERIAL); bmyz.resize(NUMBER_OF_MATERIAL);
	bmzx.resize(NUMBER_OF_MATERIAL); bmzy.resize(NUMBER_OF_MATERIAL);

	//6ñ ï™ÇÃPMLåvéZópïœêî
	//0->x0, 1->x1, 2->y0, 3->y1, 4->z0, 5->z1
	pml.resize(6);

	cout << "Matrix Resized." << endl;
}

void FDTD::InitMatrix() {
	int i, j, k, l;

	for (i = 0; i <= NX; i++) {
		for (j = 0; j <= NY; j++) {
			for (k = 0; k <= NZ; k++) {
				Ex[i][j][k] = 0.0;
				Ey[i][j][k] = 0.0;
				Ez[i][j][k] = 0.0;

				Hx[i][j][k] = 0.0;
				Hy[i][j][k] = 0.0;
				Hz[i][j][k] = 0.0;

				IDE[i][j][k] = 0;
			}
		}
	}

	cout << "Matrix Initialized." << endl;
}

void FDTD::Modeling() {
	int i, j, k;

	//IDEx,IDEy,IDEzÇÕ0Ç≈èâä˙âªÇµÇƒÇ¢ÇÈÇÃÇ≈ê^ãÛà»äOÇÃóÃàÊÇÉÇÉfÉäÉìÉOÇ∑ÇÈ
	cout << "Setting Materials." << endl;

	//ê^ãÛíÜ-0
	eps[0] = 1.0;
	sgmE[0] = 0.0;
	mu[0] = 1.0;
	sgmM[0] = 0.0;

	//êÖ-1
	eps[1] = WATER_EPS;
	sgmE[1] = WATER_SIGMA_E;
	mu[1] = WATER_MU;
	sgmM[1] = WATER_SIGMA_M;

	
	
	//ê^ãÛíÜ ID=0
	   

	//êÖ		
	for (int i = 0 + L_PML; i <= NXX + L_PML; i++) {
		for (int j = NY/2; j <= NY/2+20; j++) {
			for (int k = NZ/2-30; k <= NZ/2+30; k++) {
				IDE[i][j][k] = 1;
			}
		}
	}


	cout << "Modeling Finished." << endl;
}

void FDTD::CalcCoefficient() {
	int i, j, k;
	
	CalcCE(aex, bexy, bexz, DY, DZ);
	CalcCE(aey, beyz, beyx, DZ, DX);
	CalcCE(aez, bezx, bezy, DX, DY);
	CalcCM(amx, bmxy, bmxz, DY, DZ);
	CalcCM(amy, bmyz, bmyx, DZ, DX);
	CalcCM(amz, bmzx, bmzy, DX, DY);
	cout << "CECM Calculated." << endl;

	CalcPMLCECM(&pml[0], 0, L_PML, 0, NY, 0, NZ);
	CalcPMLCECM(&pml[1], NX - L_PML, NX, 0, NY, 0, NZ);
	CalcPMLCECM(&pml[2], 0, NX, 0, L_PML, 0, NZ);
	CalcPMLCECM(&pml[3], 0, NX, NY - L_PML, NY, 0, NZ);
	CalcPMLCECM(&pml[4], 0, NX, 0, NY, 0, L_PML);
	CalcPMLCECM(&pml[5], 0, NX, 0, NY, NZ - L_PML, NZ);
	cout << "CECM for PML Calculated." << endl;
}

void FDTD::CalcCE(vector<double> &ae, vector<double> &be1, vector<double> &be2, double d1, double d2) {
	double a;
	for (int i = 0; i < NUMBER_OF_MATERIAL; i++) {
		a = sgmE[i] * DT / (2 * eps[i] * EPSILON0);
		ae[i] = (1.0 - a) / (1.0 + a);
		be1[i] = DT / (eps[i] * EPSILON0*(1.0 + a)*d1);
		be2[i] = DT / (eps[i] * EPSILON0*(1.0 + a)*d2);
	}
}

void FDTD::CalcCM(vector<double> &am, vector<double> &bm1, vector<double> &bm2, double d1, double d2) {
	double a;
	for (int i = 0; i < NUMBER_OF_MATERIAL; i++) {
		a = sgmM[i] * DT / (2 * mu[i] * MU0);
		am[i] = (1.0 - a) / (1.0 + a);
		bm1[i] = DT / (mu[i] * MU0*(1.0 + a)*d1);
		bm2[i] = DT / (mu[i] * MU0*(1.0 + a)*d2);
	}
}

void FDTD::CalcPMLCECM(PML *pml, int nx0, int nx1, int ny0, int ny1, int nz0, int nz1) {
	int i, j, k;
	double copml = -1.528063e-4;

	//PMLÉNÉâÉXÇÃèâä˙âªÅAóÃàÊämï€
	pml->nx0 = nx0;
	pml->nx1 = nx1;
	pml->ny0 = ny0;
	pml->ny1 = ny1;
	pml->nz0 = nz0;
	pml->nz1 = nz1;
	pml->Exy.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1))); pml->Exz.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1)));
	pml->Eyz.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1))); pml->Eyx.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1)));
	pml->Ezx.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1))); pml->Ezy.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1)));
	pml->Hxy.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1))); pml->Hxz.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1)));
	pml->Hyz.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1))); pml->Hyx.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1)));
	pml->Hzx.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1))); pml->Hzy.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1)));
	pml->aexy.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1))); pml->aexz.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1)));
	pml->aeyx.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1))); pml->aeyz.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1)));
	pml->aezx.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1))); pml->aezy.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1)));
	pml->bexy.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1))); pml->bexz.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1)));
	pml->beyx.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1))); pml->beyz.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1)));
	pml->bezx.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1))); pml->bezy.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1)));
	pml->amxy.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1))); pml->amxz.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1)));
	pml->amyx.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1))); pml->amyz.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1)));
	pml->amzx.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1))); pml->amzy.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1)));
	pml->bmxy.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1))); pml->bmxz.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1)));
	pml->bmyx.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1))); pml->bmyz.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1)));
	pml->bmzx.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1))); pml->bmzy.resize(nx1 - nx0 + 1, vector<vector<double>>(ny1 - ny0 + 1, vector<double>(nz1 - nz0 + 1)));

	for (i = 0; i <= nx1 - nx0; i++) {
		for (j = 0; j <= ny1 - ny0; j++) {
			for (k = 0; k <= nz1 - nz0; k++) {
				pml->Exy[i][j][k] = pml->Exz[i][j][k] = 0;
				pml->Eyx[i][j][k] = pml->Eyz[i][j][k] = 0;
				pml->Ezx[i][j][k] = pml->Ezy[i][j][k] = 0;
				pml->Hxy[i][j][k] = pml->Hxz[i][j][k] = 0;
				pml->Hyx[i][j][k] = pml->Hyz[i][j][k] = 0;
				pml->Hzx[i][j][k] = pml->Hzy[i][j][k] = 0;
			}
		}
	}


	//åWêîåvéZ
	double smax0x = copml * REF0*(M + 1) / (L_PML*DX);
	double smax0y = copml * REF0*(M + 1) / (L_PML*DY);
	double smax0z = copml * REF0*(M + 1) / (L_PML*DZ);

	double sgmxm, sgmxe;
	double sgmym, sgmye;
	double sgmzm, sgmze;

#pragma omp parallel for
	for (k = 0; k <= nz1 - nz0 - 1; k++) {
		for (j = 0; j <= ny1 - ny0 - 1; j++) {
			for (i = 0; i <= nx1 - nx0 - 1; i++) {
				//sigmaX max ÇÃåvéZ
				if (i + nx0 < L_PML) {
					sgmxm = pow((L_PML - (i + nx0) - 0.5) / L_PML, M)*smax0x;
					sgmxe = pow(((double)L_PML - (i + nx0)) / L_PML, M)*smax0x;
				} else if (i + nx0 >= NX - L_PML) {
					sgmxm = pow(((i + nx0) - NX + L_PML + 0.5) / L_PML, M)*smax0x;
					sgmxe = pow(((double)(i + nx0) - NX + L_PML) / L_PML, M)*smax0x;
				} else {
					sgmxm = sgmxe = 0;
				}

				//sigmaY max ÇÃåvéZ
				if (j + ny0 < L_PML) {
					sgmym = pow((L_PML - (j + ny0) - 0.5) / L_PML, M)*smax0y;
					sgmye = pow(((double)L_PML - (j + ny0)) / L_PML, M)*smax0y;
				} else if (j + ny0 >= NY - L_PML) {
					sgmym = pow(((j + ny0) - NY + L_PML + 0.5) / L_PML, M)*smax0y;
					sgmye = pow(((double)(j + ny0) - NY + L_PML) / L_PML, M)*smax0y;
				} else {
					sgmym = sgmye = 0;
				}

				//sigmaZ max ÇÃåvéZ
				if (k + nz0 < L_PML) {
					sgmzm = pow((L_PML - (k + nz0) - 0.5) / L_PML, M)*smax0z;
					sgmze = pow(((double)L_PML - (k + nz0)) / L_PML, M)*smax0z;
				} else if (k + nz0 >= NZ - L_PML) {
					sgmzm = pow(((k + nz0) - NZ + L_PML + 0.5) / L_PML, M)*smax0z;
					sgmze = pow(((double)(k + nz0) - NZ + L_PML) / L_PML, M)*smax0z;
				} else {
					sgmzm = sgmze = 0;
				}

				double a;
				//ExÇÃåWêî
				a = sgmye * DT / (2 * EPSILON0);
				pml->aexy[i][j][k] = (1.0 - a) / (1.0 + a);
				pml->bexy[i][j][k] = DT / (EPSILON0*(1.0 + a)*DY);
				a = sgmze * DT / (2 * EPSILON0);
				pml->aexz[i][j][k] = (1.0 - a) / (1.0 + a);
				pml->bexz[i][j][k] = DT / (EPSILON0*(1.0 + a)*DZ);

				//EyÇÃåWêî
				a = sgmze * DT / (2 * EPSILON0);
				pml->aeyz[i][j][k] = (1.0 - a) / (1.0 + a);
				pml->beyz[i][j][k] = DT / (EPSILON0*(1.0 + a)*DZ);
				a = sgmxe * DT / (2 * EPSILON0);
				pml->aeyx[i][j][k] = (1.0 - a) / (1.0 + a);
				pml->beyx[i][j][k] = DT / (EPSILON0*(1.0 + a)*DX);

				//EzÇÃåWêî
				a = sgmxe * DT / (2 * EPSILON0);
				pml->aezx[i][j][k] = (1.0 - a) / (1.0 + a);
				pml->bezx[i][j][k] = DT / (EPSILON0*(1.0 + a)*DX);
				a = sgmye * DT / (2 * EPSILON0);
				pml->aezy[i][j][k] = (1.0 - a) / (1.0 + a);
				pml->bezy[i][j][k] = DT / (EPSILON0*(1.0 + a)*DY);

				//HxÇÃåWêî
				a = sgmym * DT / (2 * EPSILON0);
				pml->amxy[i][j][k] = (1.0 - a) / (1.0 + a);
				pml->bmxy[i][j][k] = DT / (MU0*(1.0 + a)*DY);
				a = sgmzm * DT / (2 * EPSILON0);
				pml->amxz[i][j][k] = (1.0 - a) / (1.0 + a);
				pml->bmxz[i][j][k] = DT / (MU0*(1.0 + a)*DZ);

				//HyÇÃåWêî
				a = sgmxm * DT / (2 * EPSILON0);
				pml->amyx[i][j][k] = (1.0 - a) / (1.0 + a);
				pml->bmyx[i][j][k] = DT / (MU0*(1.0 + a)*DX);
				a = sgmzm * DT / (2 * EPSILON0);
				pml->amyz[i][j][k] = (1.0 - a) / (1.0 + a);
				pml->bmyz[i][j][k] = DT / (MU0*(1.0 + a)*DZ);

				//HzÇÃåWêî
				a = sgmym * DT / (2 * EPSILON0);
				pml->amzy[i][j][k] = (1.0 - a) / (1.0 + a);
				pml->bmzy[i][j][k] = DT / (MU0*(1.0 + a)*DY);
				a = sgmxm * DT / (2 * EPSILON0);
				pml->amzx[i][j][k] = (1.0 - a) / (1.0 + a);
				pml->bmzx[i][j][k] = DT / (MU0*(1.0 + a)*DX);
			}
		}
	}
}

void FDTD::SourceE(double t) {
	int j = NY/2-20;
	for (int i = 0; i <= NX; i++) {
		for (int k = 0; k <= NZ; k++) {
			Ez[i][j][k] = E_WAVE_AMPLITUDE * sin(2 * PI*WAVE_FREQUENCY*t);
			//Ey[i][j][k] = E_WAVE_AMPLITUDE * sin(2 * PI*WAVE_FREQUENCY*t)*sin(PI / 6);
			//Ez[i][j][k] = E_WAVE_AMPLITUDE * sin(2 * PI*WAVE_FREQUENCY*t)*cos(PI / 6);
		}
	}
}

void FDTD::SourceH(double t) {
	int j = NY / 2 - 20;
	for (int i = 0; i <= NX; i++) {
		for (int k = 0; k <= NZ; k++) {
			Hx[i][j][k] = H_WAVE_AMPLITUDE * sin(2 * PI*WAVE_FREQUENCY*t);
			//Hx[i][j][k] = H_WAVE_AMPLITUDE * sin(2 * PI*WAVE_FREQUENCY*t)*sin(PI/6);
		}
	}
}

void FDTD::CalcEField() {
	int i, j, k;
	int id;

	//Ex
#pragma omp parallel for
	for (k = 1; k <= NZ - 1; k++) {
		for (j = 1; j <= NY - 1; j++) {
			for (i = 0; i <= NX - 1; i++) {
				id = IDE[i][j][k];
				Ex[i][j][k] = aex[id] * Ex[i][j][k] + bexy[id] * (Hz[i][j][k] - Hz[i][j - 1][k]) - bexz[id] * (Hy[i][j][k] - Hy[i][j][k - 1]);
			}
		}
	}

	//Ey
#pragma omp parallel for
	for (k = 1; k <= NZ - 1; k++) {
		for (j = 0; j <= NY - 1; j++) {
			for (i = 1; i <= NX - 1; i++) {
				id = IDE[i][j][k];
				Ey[i][j][k] = aey[id] * Ey[i][j][k] + beyz[id] * (Hx[i][j][k] - Hx[i][j][k - 1]) - beyx[id] * (Hz[i][j][k] - Hz[i - 1][j][k]);
			}
		}
	}

	//Ez
#pragma omp parallel for
	for (k = 0; k <= NZ - 1; k++) {
		for (j = 1; j <= NY - 1; j++) {
			for (i = 1; i <= NX - 1; i++) {
				id = IDE[i][j][k];
				Ez[i][j][k] = aez[id] * Ez[i][j][k] + bezx[id] * (Hy[i][j][k] - Hy[i - 1][j][k]) - bezy[id] * (Hx[i][j][k] - Hx[i][j - 1][k]);
			}
		}
	}

	//cout << "Updated E Field." << endl;
}

void FDTD::CalcHField() {
	int i, j, k;
	int id;

	//Hx
#pragma omp parallel for
	for (k = 0; k <= NZ - 1; k++) {
		for (j = 0; j <= NY - 1; j++) {
			for (i = 1; i <= NX - 1; i++) {
				id = IDE[i][j][k];
				Hx[i][j][k] = amx[id] * Hx[i][j][k] - bmxy[id] * (Ez[i][j + 1][k] - Ez[i][j][k]) + bmxz[id] * (Ey[i][j][k + 1] - Ey[i][j][k]);
			}
		}
	}

	//Hy
#pragma omp parallel for
	for (k = 0; k <= NZ - 1; k++) {
		for (j = 1; j <= NY - 1; j++) {
			for (i = 0; i <= NX - 1; i++) {
				id = IDE[i][j][k];
				Hy[i][j][k] = amy[id] * Hy[i][j][k] - bmyz[id] * (Ex[i][j][k + 1] - Ex[i][j][k]) + bmyx[id] * (Ez[i + 1][j][k] - Ez[i][j][k]);
			}
		}
	}

	//Hz
#pragma omp parallel for
	for (k = 1; k <= NZ - 1; k++) {
		for (j = 0; j <= NY - 1; j++) {
			for (i = 0; i <= NX - 1; i++) {
				id = IDE[i][j][k];
				Hz[i][j][k] = amz[id] * Hz[i][j][k] - bmzx[id] * (Ey[i + 1][j][k] - Ey[i][j][k]) + bmzy[id] * (Ex[i][j + 1][k] - Ex[i][j][k]);
			}
		}
	}
	//cout << "Updated H Field." << endl;
}

void FDTD::CalcPMLEField(PML *pml) {
	int i, j, k;

	//Ex
#pragma omp parallel for
	for (k = 1; k <= pml->nz1 - pml->nz0 - 1; k++) {
		for (j = 1; j <= pml->ny1 - pml->ny0 - 1; j++) {
			for (i = 0; i <= pml->nx1 - pml->nx0 - 1; i++) {
				pml->Exy[i][j][k] = pml->aexy[i][j][k] * pml->Exy[i][j][k]
					+ pml->bexy[i][j][k] * (Hz[i + pml->nx0][j + pml->ny0][k + pml->nz0] - Hz[i + pml->nx0][j - 1 + pml->ny0][k + pml->nz0]);
				pml->Exz[i][j][k] = pml->aexz[i][j][k] * pml->Exz[i][j][k]
					+ pml->bexz[i][j][k] * (Hy[i + pml->nx0][j + pml->ny0][k - 1 + pml->nz0] - Hy[i + pml->nx0][j + pml->ny0][k + pml->nz0]);
				Ex[i + pml->nx0][j + pml->ny0][k + pml->nz0] = pml->Exy[i][j][k] + pml->Exz[i][j][k];
			}
		}
	}

	//Ey
#pragma omp parallel for
	for (k = 1; k <= pml->nz1 - pml->nz0 - 1; k++) {
		for (j = 0; j <= pml->ny1 - pml->ny0 - 1; j++) {
			for (i = 1; i <= pml->nx1 - pml->nx0 - 1; i++) {
				pml->Eyz[i][j][k] = pml->aeyz[i][j][k] * pml->Eyz[i][j][k]
					+ pml->beyz[i][j][k] * (Hx[i + pml->nx0][j + pml->ny0][k + pml->nz0] - Hx[i + pml->nx0][j + pml->ny0][k - 1 + pml->nz0]);
				pml->Eyx[i][j][k] = pml->aeyx[i][j][k] * pml->Eyx[i][j][k]
					+ pml->beyx[i][j][k] * (Hz[i - 1 + pml->nx0][j + pml->ny0][k + pml->nz0] - Hz[i + pml->nx0][j + pml->ny0][k + pml->nz0]);
				Ey[i + pml->nx0][j + pml->ny0][k + pml->nz0] = pml->Eyz[i][j][k] + pml->Eyx[i][j][k];
			}
		}
	}

	//Ez
#pragma omp parallel for
	for (k = 0; k <= pml->nz1 - pml->nz0 - 1; k++) {
		for (j = 1; j <= pml->ny1 - pml->ny0 - 1; j++) {
			for (i = 1; i <= pml->nx1 - pml->nx0 - 1; i++) {
				pml->Ezx[i][j][k] = pml->aezx[i][j][k] * pml->Ezx[i][j][k]
					+ pml->bezx[i][j][k] * (Hy[i + pml->nx0][j + pml->ny0][k + pml->nz0] - Hy[i - 1 + pml->nx0][j + pml->ny0][k + pml->nz0]);
				pml->Ezy[i][j][k] = pml->aezy[i][j][k] * pml->Ezy[i][j][k]
					+ pml->bezy[i][j][k] * (Hx[i + pml->nx0][j - 1 + pml->ny0][k + pml->nz0] - Hx[i + pml->nx0][j + pml->ny0][k + pml->nz0]);
				Ez[i + pml->nx0][j + pml->ny0][k + pml->nz0] = pml->Ezx[i][j][k] + pml->Ezy[i][j][k];
			}
		}
	}

	//cout << "Updated PML E Field." << endl;
}

void FDTD::CalcPMLHField(PML *pml) {
	int i, j, k;

	//Hx
#pragma omp parallel for
	for (k = 0; k <= pml->nz1 - pml->nz0 - 1; k++) {
		for (j = 0; j <= pml->ny1 - pml->ny0 - 1; j++) {
			for (i = 1; i <= pml->nx1 - pml->nx0 - 1; i++) {
				pml->Hxy[i][j][k] = pml->amxy[i][j][k] * pml->Hxy[i][j][k]
					+ pml->bmxy[i][j][k] * (Ez[i + pml->nx0][j + pml->ny0][k + pml->nz0] - Ez[i + pml->nx0][j + 1 + pml->ny0][k + pml->nz0]);
				pml->Hxz[i][j][k] = pml->amxz[i][j][k] * pml->Hxz[i][j][k]
					+ pml->bmxz[i][j][k] * (Ey[i + pml->nx0][j + pml->ny0][k + 1 + pml->nz0] - Ey[i + pml->nx0][j + pml->ny0][k + pml->nz0]);
				Hx[i + pml->nx0][j + pml->ny0][k + pml->nz0] = pml->Hxy[i][j][k] + pml->Hxz[i][j][k];
			}
		}
	}

	//Hy
#pragma omp parallel for
	for (k = 0; k <= pml->nz1 - pml->nz0 - 1; k++) {
		for (j = 1; j <= pml->ny1 - pml->ny0 - 1; j++) {
			for (i = 0; i <= pml->nx1 - pml->nx0 - 1; i++) {
				pml->Hyz[i][j][k] = pml->amyz[i][j][k] * pml->Hyz[i][j][k]
					+ pml->bmyz[i][j][k] * (Ex[i + pml->nx0][j + pml->ny0][k + pml->nz0] - Ex[i + pml->nx0][j + pml->ny0][k + 1 + pml->nz0]);
				pml->Hyx[i][j][k] = pml->amyx[i][j][k] * pml->Hyx[i][j][k]
					+ pml->bmyx[i][j][k] * (Ez[i + 1 + pml->nx0][j + pml->ny0][k + pml->nz0] - Ez[i + pml->nx0][j + pml->ny0][k + pml->nz0]);
				Hy[i + pml->nx0][j + pml->ny0][k + pml->nz0] = pml->Hyz[i][j][k] + pml->Hyx[i][j][k];
			}
		}
	}

	//Hz
#pragma omp parallel for
	for (k = 1; k <= pml->nz1 - pml->nz0 - 1; k++) {
		for (j = 0; j <= pml->ny1 - pml->ny0 - 1; j++) {
			for (i = 0; i <= pml->nx1 - pml->nx0 - 1; i++) {
				pml->Hzx[i][j][k] = pml->amzx[i][j][k] * pml->Hzx[i][j][k]
					+ pml->bmzx[i][j][k] * (Ey[i + pml->nx0][j + pml->ny0][k + pml->nz0] - Ey[i + 1 + pml->nx0][j + pml->ny0][k + pml->nz0]);
				pml->Hzy[i][j][k] = pml->aezy[i][j][k] * pml->Hzy[i][j][k]
					+ pml->bmzy[i][j][k] * (Ex[i + pml->nx0][j + 1 + pml->ny0][k + pml->nz0] - Ex[i + pml->nx0][j + pml->ny0][k + pml->nz0]);
				Hz[i + pml->nx0][j + pml->ny0][k + pml->nz0] = pml->Hzx[i][j][k] + pml->Hzy[i][j][k];
			}
		}
	}


	//cout << "Updated PML H Field." << endl;
}

void FDTD::StartRepetition() {
	int i, j, k, n;
	double t = DT;
	char file_name[20];
	Output* EOutput1;
	Output *EOutput2;
	cout << "Repetition Start." << endl;

	for (n = 0; n <= NT; n++) {
		CalcEField();
		for (i = 0; i < 6; i++) {
			CalcPMLEField(&pml[i]);
		}
		SourceE(t);
		t = t + DT / 2.0;
		CalcHField();
		for (i = 0; i < 6; i++) {
			CalcPMLHField(&pml[i]);
		}
		SourceH(t);
		t = t + DT / 2.0;


		cout << "n=" << setw(4) << setfill('0') << n;
		cout << " (t=" << scientific << setprecision(2) << t << ")" << endl;

		if (n % 50 == 0) {
			sprintf(file_name, "E-n=%03d,t=%4.2e.csv", n, t);
			//EOutput1 = new Output(file_name);
			sprintf(file_name, "E-n=%03d,t=%4.2e.csv", n, t);
			EOutput2 = new Output(file_name);

			for (i = L_PML; i <= NX - L_PML; i++) {
				for (j = L_PML; j <= NY - L_PML; j++) {
					for (k = L_PML; k <= NZ - L_PML; k++) {
						//EOutput1->Add(i - L_PML, j - L_PML, k - L_PML, Ez[i][j][k]);
						EOutput2->Add((i - L_PML)*DX*1e-6, (j - L_PML)*DY*1e-6, (k - L_PML)*DZ*1e-6, Ez[i][j][k]);
					}
					//EOutput1->Add();
					//EOutput2->Add();
				}
			}
			//delete EOutput1;
			delete EOutput2;
		}

	}
}




