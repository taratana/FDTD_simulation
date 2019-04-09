#include "FDTD.h"
#include<fstream>
#include <ctime>
using namespace std;

class StopWatch {
private:
	time_t  st_t, fin_t;  /*時間計測用（実時間） */
	clock_t st_ct, fin_ct; /*時間計測用（CPU時間）*/

public:
	void start() { /*時間計測開始*/
		st_ct = clock();
		time(&st_t);
	}
	void stop() { /*時間計測終了*/
		fin_ct = clock();
		time(&fin_t);
	}
	double getTime() {
		return difftime(fin_t, st_t);
	}
};

int main(void) {
	StopWatch record_time;
	record_time.start();
	ofstream of_log("output.txt");

	FDTD fdtd;

	fdtd.DeterminMatrixSize();

	fdtd.InitMatrix();

	fdtd.Modeling();
	
	fdtd.CalcCoefficient();

	fdtd.StartRepetition();

	record_time.stop();
	double elapsed = record_time.getTime();
	
	of_log << "<Number of Partitions excluding PML>-(NXX,NYY,NZZ)" << endl << "(" << NXX << "," << NYY << "," << NZZ << ")" << "\n\n";
	of_log << "<PML Number>-L" << endl << L_PML << "\n\n";
	of_log << "<DT>" << endl << DT << "\n\n";
	of_log << "<Time Step>-NT " << endl << NT << "\n\n";
	of_log << "<Cell Size>-(DX,DY,DZ)"<<endl<< "(" << DX << "," << DY << "," << DZ <<")"<< "\n\n";
	of_log << "<Material Constants>" << endl << "water-(εr,μr,σe,σm,)" << "\n";
	of_log << "(" << WATER_EPS << "," << WATER_MU << "," << WATER_SIGMA_E << "," << WATER_SIGMA_M << ")" << "\n\n";
	of_log << "<Time Elapased>" << endl << (int)elapsed / 3600 << "hours " << (int)elapsed % 3600 / 60.0 << "min" << "\n\n";
}
