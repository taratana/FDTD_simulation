#include "FDTD.h"
#include<fstream>
#include <ctime>
using namespace std;

class StopWatch {
private:
	time_t  st_t, fin_t;  /*ΤvͺpiΐΤj */
	clock_t st_ct, fin_ct; /*ΤvͺpiCPUΤj*/

public:
	void start() { /*ΤvͺJn*/
		st_ct = clock();
		time(&st_t);
	}
	void stop() { /*ΤvͺIΉ*/
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
	of_log << "<Material Constants>" << endl << "water-(Γr,Κr,Πe,Πm,)" << "\n";
	of_log << "(" << WATER_EPS << "," << WATER_MU << "," << WATER_SIGMA_E << "," << WATER_SIGMA_M << ")" << "\n\n";
	of_log << "<Time Elapased>" << endl << (int)elapsed / 3600 << "hours " << (int)elapsed % 3600 / 60.0 << "min" << "\n\n";
}
