#include "GraphR.h"


GraphR::GraphR(string name){
	fp.open(name.c_str(),ios::out);
	if(fp.fail()){
		std::cout<<"file open error\n";
	}
	typeflag = 0;
}






void GraphR::Add(double x,double y,double z,double V){
	switch(typeflag){
		case 0:	//初回書き込み
			typeflag = 1;
			fp<<"データ形式,3\n";
			fp<<"memo\n";
			fp<<"x,y,z,V\n";
			fp<<x<<","<<y<<","<<z<<","<<V<<"\n";
			break;
		case 1:
			fp<<x<<","<<y<<","<<z<<","<<V<<"\n";
			break;
		default:
			fp<<"you have some mistekes"<<"\n";
			break;
	}
}


void GraphR::Add(double x,double y,double z,double Vx,double Vy,double Vz){
	switch(typeflag){
		case 0:	//初回書き込み
			typeflag = 2;
			fp<<"データ形式,5\n";
			fp<<"memo\n";
			fp<<"x,y,z,Vx,Vy,Vz\n";
			fp<<x<<","<<y<<","<<z<<","<<Vx<<","<<Vy<<","<<Vz<<"\n";
			break;
		case 2:
			fp<<x<<","<<y<<","<<z<<","<<Vx<<","<<Vy<<","<<Vz<<"\n";
			break;
		default:
			fp<<"you have some mistekes"<<"\n";
			break;	
	}
}


void GraphR::Add(double x, double y, double V) {
	switch (typeflag) {
	case 0:	//初回書き込み
		typeflag = 3;
		fp << "データ形式,2\n";
		fp << "memo\n";
		fp << "x,y,V\n";
		fp << x << "," << y << "," << V << "\n";
		break;
	case 3:
		fp << x << "," << y << "," << V << "\n";
		break;
	default:
		fp << "you have some mistekes" << "\n";
		break;
	}
}

void GraphR::close(){
	typeflag = 0;
	fp.close();
}