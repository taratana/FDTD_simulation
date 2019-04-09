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
	//コンストラクタ
	GraphR(string name);

	//メソッド
	void Add(double x,double y,double z,double V);
	void Add(double x,double y,double z,double Vx,double Vy,double Vz);
	void Add(double x, double y, double V);
	void close();
};