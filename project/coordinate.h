#pragma once
#include<iostream>
#include<cmath>
#define PI 3.141592653589793
#define rou (180/3.141592653589793)

using namespace std;

// 地球椭球参数。
struct EARTH
{
	double a = 0;
	double f = 0;
}const WGS84 = { 6378137.0 ,1.0 / 298.257223563 }, CGCS2K = { 6378137.0 ,1.0 / 298.257222101 };

struct XYZ
{
	double x;
	double y;
	double z;
	XYZ(double X, double Y, double Z)
	{
		x = X; y = Y; z = Z;
	}
	XYZ()
	{
		x = 0; y = 0; z = 0;
	};
	XYZ operator-(const XYZ& b) const;
	XYZ operator+(const XYZ& b) const;
};

struct BLH
{
	double b = 0;
	double l = 0;
	double h = 0;
	void BLdeg()// 弧度转角度
	{
		this->b = this->b * 180 / PI;
		this->l = this->l * 180 / PI;
	}
	void BLrad()// 角度转弧度
	{
		this->b = this->b / 180 * PI;
		this->l = this->l / 180 * PI;
	}
	BLH(double bb, double ll, double hh)
	{
		b = bb; l = ll; h = hh;
	}
	BLH() { b = 0; l = 0; h = 0; }
};

struct NEU
{
	double n = 0;
	double e = 0;
	double u = 0;
};

void BLH2XYZ(EARTH Ea, BLH* blh, XYZ* xyz);
int XYZ2BLH(EARTH Ea, BLH* blh, XYZ* xyz);
int XYZ2NEU(EARTH Ea, XYZ* xyz, XYZ* xyz0, NEU* neu);
double Dist(XYZ* xyz1, XYZ* xyz2);
double VerticalAng(XYZ* xyz1, XYZ* xyz0);
XYZ XYZAdd(XYZ a, XYZ b);
