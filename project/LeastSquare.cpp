#include"decoding.h"
#include"Atomsphere.h"
using namespace std;

/*******************************
计算距离的函数，根据解算得到的卫星位置和速度
该函数对GPSInfo和BDSInfo都进行了重载，也就是都可以计算其中的卫星到测站的距离
第二个参数CMatrix X是坐标的列向量，只要前三个元素是坐标即可
***********************************/
double CalDist(GPSInfo e, CMatrix X)
{
	if (X.length < 3)
	{
		cerr << "Distance Error: Incorrect length" << endl;
		return 0;
	}
	double D;
	D = sqrt((e.X - X.value[0]) * (e.X - X.value[0]) +
		(e.Y - X.value[1]) * (e.Y - X.value[1]) + (e.Z - X.value[2]) * (e.Z - X.value[2]));
	return D;
}
double CalDist(BDSInfo e, CMatrix X)
{
	if (X.length < 3)
	{
		cerr << "Distance Error: Incorrect length" << endl;
		return 0;
	}
	double D;
	D = sqrt((e.X - X.value[0]) * (e.X - X.value[0]) +
		(e.Y - X.value[1]) * (e.Y - X.value[1]) + (e.Z - X.value[2]) * (e.Z - X.value[2]));
	return D;
}

double CalDist(XYZ p, XYZ q)
{
	double D = sqrt((p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y) + (p.z - q.z) * (p.z - q.z));
	return D;
}
/*****************************
计算前三个分量的模长
即计算XYZ的模长
****************************/
double CalNorm3(CMatrix C)
{
	double* D = C.value;
	double R = sqrt(D[0] * D[0] + D[1] * D[1] + D[2] * D[2]);
	return R;
}

/*
最小二乘迭代函数，输入的信息有：
	单历元结构体
	空的B,W,X0,dX矩阵，该函数会进行最小二乘迭代
最终这几个矩阵都会自行更新到迭代收敛结果，如果出现问题，将会有提示且返回0值
*/
int LeastSquareIteration(ONEpoch* e, CMatrix& B, CMatrix& W, CMatrix& XR0, CMatrix& d_X)
{
	int num = 0;
	while (CalNorm3(d_X) > 1e-5)
	{
		num++;
		if (num > 1000)
		{
			cout << "警告：最小二乘迭代陷入死循环" << endl;
			return 0;
		}
		XYZ xyz0(XR0.value[0], XR0.value[1], XR0.value[2]);
		CMatrix Nbb(5, 5);
		for (int i = 1; i < e->nGvalid + 1; i++)
		{
			GPSInfo G1 = e->gInfo[i];
			XYZ xyz1(G1.X, G1.Y, G1.Z);
			NEU neu1;
			XYZ2NEU(WGS84, &xyz1, &xyz0, &neu1);
			double Trop = HopeField(St, neu1);//计算对流层改正参数（霍普菲尔德模型）
			e->gInfo[i].Trop = Trop;
			double L1 = e->gObs[i].raw[0].psr; //- Trop; //+C_Light * e->gInfo[i].dt - XR0.value[3];
			double L2 = e->gObs[i].raw[1].psr;// -Trop; //+C_Light * e->gInfo[i].dt - XR0.value[3];

			//为系数矩阵赋值
			double distG = CalDist(e->gInfo[i], XR0);
			B.value[5 * (i - 1)] = (xyz0.x - G1.X) / distG;
			B.value[5 * (i - 1) + 1] = (xyz0.y - G1.Y) / distG;
			B.value[5 * (i - 1) + 2] = (xyz0.z - G1.Z) / distG;
			B.value[5 * (i - 1) + 3] = 1;
			B.value[5 * (i - 1) + 4] = 0;
			double Ig = GIonFree(L1, L2);//无电离层组合观测值
			W.value[i - 1] = Ig - distG - Trop + C_Light * e->gInfo[i].dt - XR0.value[3];
		}
		for (int i = 1; i < e->nBvalid + 1; i++)
		{
			int btemp = i + e->nGvalid;
			BDSInfo B1 = e->bInfo[i];
			XYZ xyz1(B1.X, B1.Y, B1.Z);
			NEU neu1;
			XYZ2NEU(CGCS2K, &xyz1, &xyz0, &neu1);						//计算站心坐标系
			double Trop = HopeField(St, neu1);										//计算对流层改正参数（霍普菲尔德模型）
			e->bInfo[i].Trop = Trop;
			double L1 = e->bObs[i].raw[0].psr - e->bInfo[i].tgd1 * C_Light;// -Trop + C_Light * e->bInfo[i].dt - XR0.value[4];
			double L3 = e->bObs[i].raw[1].psr;// -Trop + C_Light * e->bInfo[i].dt - XR0.value[4];
			double distB = CalDist(e->bInfo[i], XR0);
			B.value[5 * (btemp - 1)] = (xyz0.x - B1.X) / distB;
			B.value[5 * (btemp - 1) + 1] = (xyz0.y - B1.Y) / distB;
			B.value[5 * (btemp - 1) + 2] = (xyz0.z - B1.Z) / distB;
			B.value[5 * (btemp - 1) + 3] = 0;
			B.value[5 * (btemp - 1) + 4] = 1;
			double Ib = BIonFree(L1, L3);
			W.value[btemp - 1] = Ib - distB - Trop + C_Light * e->bInfo[i].dt - XR0.value[4];
		}


		Nbb = (B.transpose() * B).inv();
		d_X = Nbb * B.transpose() * W;
		XR0 = XR0 + d_X;
		//MatrixShow(Nbb);

	}
	return 1;
}