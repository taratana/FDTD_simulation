#include "Output.h"


Output::Output(string name){
	fp.open(name.c_str(),ios::out);
	if(fp.fail()){
		std::cout<<"file open error\n";
	}
	typeflag = 0;
}






void Output::Add(double x,double y,double z,double V){
	switch(typeflag){
		case 0:	//���񏑂�����
			typeflag = 1;
			fp<<"#gnuplot�œǂݍ��ނ��ƂŃv���b�g�ł��܂��B\n";
			fp << "#���̂Ƃ���every���g���ĊJ�n�u���b�N�̐錾�Ȃǂ��s���Ƃ��܂��X�N���v�g���ł��܂��B\n";
			fp<<"#x,y,z,V\n";
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


void Output::Add(double x,double y,double z,double Vx,double Vy,double Vz){
	switch(typeflag){
		case 0:	//���񏑂�����
			typeflag = 2;
			fp<<"�f�[�^�`��,5\n";
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


void Output::Add(double x, double y, double V) {
	switch (typeflag) {
	case 0:	//���񏑂�����
		typeflag = 3;
		fp << "�f�[�^�`��,2\n";
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

void Output::Add() {
	fp << endl;
}

void Output::close(){
	typeflag = 0;
	fp.close();
}