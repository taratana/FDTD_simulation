#pragma once
#include<fstream>
#include<iostream>
#include<string>

using namespace std;

class GraphR{
private:
	ofstream fp;
	int typeflag;
public:
	//�R���X�g���N�^
	GraphR(string name);

	//���\�b�h
	void Add(double x,double y,double z,double V);
	void Add(double x,double y,double z,double Vx,double Vy,double Vz);
	void Add(double x, double y, double V);
	void close();
};