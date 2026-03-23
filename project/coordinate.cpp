#include"coordinate.h"
/*************************************
BLHToXYZ
目的：大地坐标到笛卡尔坐标的转换算法

参数：
a   椭球长半轴
f   椭球扁率
blh 大地坐标系坐标
xyz ecef坐标系坐标
**************************************/
void BLH2XYZ(EARTH Ea, BLH* blh, XYZ* xyz)
{
	double a = Ea.a; double f = Ea.f;
	double N, e2, B, L, H;
	B = blh->b;
	L = blh->l;
	H = blh->h;

	e2 = 2 * f - f * f;
	N = a / (sqrt(1 - e2 * sin(B) * sin(B)));
	xyz->x = (N + H) * cos(B) * cos(L);
	xyz->y = (N + H) * cos(B) * sin(L);
	xyz->z = (N * (1 - e2) + H) * sin(B);
}

/*************************************
XYZ2BLH
目的：笛卡尔坐标到大地坐标的转换算法

参数：
a   椭球长半轴
f   椭球扁率
B   大地经度（单位：弧度）
blh 大地坐标系坐标
xyz ecef坐标系坐标

返回值：1=正常 0=死循环
**************************************/
int XYZ2BLH(EARTH Ea, BLH* blh, XYZ* xyz)
{
	double a = Ea.a; double f = Ea.f;
	int j;
	double delta_Z[2], N, e2, X, Y, Z;
	X = xyz->x;
	Y = xyz->y;
	Z = xyz->z;

	if ((X * X + Y * Y + Z * Z) < 1e-8)
	{
		blh->b = 0.0;
		blh->l = 0.0;
		blh->h = -a;
		return 0;
	}
	e2 = 2 * f - f * f;
	delta_Z[0] = e2 * Z;

	blh->l = atan2(Y, X);
	for (j = 0; j < 10; j++)
	{
		blh->b = atan2(Z + delta_Z[0], sqrt(X * X + Y * Y));
		N = a / (sqrt(1 - e2 * sin(blh->b) * sin(blh->b)));
		blh->h = sqrt(X * X + Y * Y + (Z + delta_Z[0]) * (Z + delta_Z[0])) - N;
		delta_Z[1] = N * e2 * (Z + delta_Z[0]) / (sqrt(X * X + Y * Y + (Z + delta_Z[0]) * (Z + delta_Z[0])));

		if (abs(delta_Z[1] - delta_Z[0]) < 1e-10) break;
		else delta_Z[0] = delta_Z[1];
	}

	if (abs(delta_Z[1] - delta_Z[0]) >= 1e-10)
	{
		cerr << "XYZ2BLH进入死循环！" << endl;
		return 0;
	}
	else
	{
		blh->b = atan2(Z + delta_Z[1], sqrt(X * X + Y * Y));
		N = a / (sqrt(1 - e2 * sin(blh->b) * sin(blh->b)));
		blh->h = sqrt(X * X + Y * Y + (Z + delta_Z[1]) * (Z + delta_Z[1])) - N;
		return 1;
	}
}
/*************************
XYZ转NEU坐标系的程序，需要的参数有
地球参数Ea
原点坐标xyz0
待测点坐标xyz
站心坐标（输出的）neu
**************************/
int XYZ2NEU(EARTH Ea, XYZ* xyz, XYZ* xyz0, NEU* neu)
{
	BLH blh;
	int R = XYZ2BLH(Ea, &blh, xyz0);
	if (R == 0)return R;
	double N = -sin(blh.b) * cos(blh.l) * (xyz->x - xyz0->x) - sin(blh.b) * sin(blh.l) * (xyz->y - xyz0->y) + cos(blh.b) * (xyz->z - xyz0->z);
	double E = -sin(blh.l) * (xyz->x - xyz0->x) + cos(blh.l) * (xyz->y - xyz0->y);
	double U = cos(blh.b) * cos(blh.l) * (xyz->x - xyz0->x) + cos(blh.b) * sin(blh.l) * (xyz->y - xyz0->y) + sin(blh.b) * (xyz->z - xyz0->z);
	neu->e = E;
	neu->u = U;
	neu->n = N;
	return 1;
}

/***************************
计算距离的函数，最后返回的是两者的距离
输入参数是两点的坐标，返回参数是距离值
************************/
double Dist(XYZ* xyz1, XYZ* xyz2)
{
	double dx = xyz1->x - xyz2->x;
	double dy = xyz1->y - xyz2->y;
	double dz = xyz1->z - xyz2->z;
	double d = sqrt(dx * dx + dy * dy + dz * dz);
	return d;
}

/************************
计算高度角的函数，注意基站位置应该放在后面
这里会直接转换站心坐标系，所以一定要保证两点的顺序
最后返回的数值是弧度制的高度角（有正有负）
*************************/
double VerticalAng(XYZ* xyz1, XYZ* xyz0)
{
	NEU neu1;
	XYZ2NEU(WGS84, xyz1, xyz0, &neu1);
	double ne = sqrt(neu1.e * neu1.e + neu1.n * neu1.n);
	double ang = atan(neu1.u / ne);
	return ang;
}

/*
	XYZ位置坐标向量相加的函数
*/
XYZ XYZAdd(XYZ a, XYZ b)
{
	XYZ Sum(a.x + b.x, a.y + b.y, a.z + b.z);
	return Sum;
}

XYZ XYZ::operator-(const XYZ& b) const
{
	XYZ Res(x - b.x, y - b.y, z - b.z);
	return Res;
}

XYZ XYZ::operator+(const XYZ& b) const
{
	XYZ Res(x + b.x, y + b.y, z + b.z);
	return Res;
}