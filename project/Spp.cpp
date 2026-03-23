#include"decoding.h"
#include"Atomsphere.h"
#include<fstream>
#include<iomanip>
#include<vector>
using namespace std;

/*************************************
CheckEpoch函数，主要功能是检查该历元的结构体能不能正常解算
不可用返回0，可用返回1
它可以用来检查星历的可用性
*************************************/
int CheckEpoch(ONEpoch* e)
{
	if (e->bEphValid == false || e->gEphValid == false)//如果星历不可用
	{
		cout << "星历不可用" << endl;
		return 0;
	}
	else if ((e->nGvalid + e->nBvalid) < 5 || e->nGvalid == 0 || e->nBvalid == 0) //如果可用的卫星总数小于5，不满足多余观测要求
	{
		cout << "观测值过少" << endl;
		return 0;
	}
	else return 1;
	return 0;
}

/*
单点定位函数：
其需要输入单个历元观测结构体e以及SPP结果结构体SPPReslut
进行最小二乘迭代解算，并输出结果
*/
int StdPos(ONEpoch* e, SPPResult* sppr)
{
	double x0[] = { -2267804.5263,5009342.3723,3220991.8632,0.089,0.123 };//为测站赋给一个大致的初值
	CMatrix XR0(x0, 5);
	GetSatInfo(e);//计算信号发出时刻的卫星位置、钟差和钟速
	XYZ xyz0(XR0.value[0], XR0.value[1], XR0.value[2]);

	//初始化局部变量

	CMatrix B(e->nBvalid + e->nGvalid, 5), W(e->nBvalid + e->nGvalid, 1);
	double a[] = { 100,100,100,100,100 };//赋予它一个初值
	CMatrix d_X(a, 5);
	//迭代进行最小二乘解算，阈值设为模长的误差矩阵小于1e-5
	LeastSquareIteration(e, B, W, XR0, d_X);
	XYZ xyzR(XR0.value[0], XR0.value[1], XR0.value[2]);
	DetectAngH(&xyzR, e, 10.0);
	B = CMatrix(e->nBvalid + e->nGvalid, 5); W = CMatrix(e->nBvalid + e->nGvalid, 1);
	d_X = CMatrix(a, 5);
	LeastSquareIteration(e, B, W, XR0, d_X);
	CMatrix V = B * d_X + (-W);
	double Sig0 = sqrt((V.transpose() * V).value[0] / (e->nBvalid + e->nGvalid - 5));
	//PsoteriorOutlier(e, V, Sig0);
	//B = CMatrix(e->nBvalid + e->nGvalid, 5); W = CMatrix(e->nBvalid + e->nGvalid, 1);
	//d_X = CMatrix(a, 5);
	LeastSquareIteration(e, B, W, XR0, d_X);
	//////对验后粗差进行剔除
	//V = B * d_X + (-W);
	//Sig0 = sqrt((V.transpose()*V).value[0] / (e->nBvalid + e->nGvalid - 5));
	CMatrix Q = (B.transpose() * B).inv();
	sppr->Sig0 = Sig0;
	sppr->Gtime = e->Gtime;
	sppr->X = XR0.value[0]; sppr->Y = XR0.value[1]; sppr->Z = XR0.value[2];
	sppr->dtG = XR0.value[3]; sppr->dtB = XR0.value[4];
	sppr->nBvalid = e->nBvalid; sppr->nGvalid = e->nGvalid;
	double* v = V.value;
	sppr->PDOP = sqrt(Q.find(1, 1) * Q.find(1, 1) + Q.find(2, 2) * Q.find(2, 2) + Q.find(3, 3) * Q.find(3, 3));
	//for (int s = 0; s < XR0.length; s++)cout << XR0.value[s]<<endl;
	XYZ xyztemp(sppr->X, sppr->Y, sppr->Z);
	BLH blh;
	XYZ2BLH(WGS84, &blh, &xyztemp);
	sppr->blh = blh;
	return 1;
}
/*********************************
计算无电离层组合的函数
GPS和BDS参数不同，采用两个函数进行解算
************************************/
double GIonFree(double a, double b)
{
	double Res;
	Res = (FG1_GPS * FG1_GPS) / (FG1_GPS * FG1_GPS - FG2_GPS * FG2_GPS) * a
		- (FG2_GPS * FG2_GPS) / (FG1_GPS * FG1_GPS - FG2_GPS * FG2_GPS) * b;
	return Res;
}

double BIonFree(double a, double b)
{
	double Res;
	Res = (FB1_BDS * FB1_BDS) / (FB1_BDS * FB1_BDS - FB3_BDS * FB3_BDS) * a
		- (FB3_BDS * FB3_BDS) / (FB1_BDS * FB1_BDS - FB3_BDS * FB3_BDS) * b;
	return Res;
}
/***************************************
单点测速程序，需要运行了定位函数后继续使用
单点测速需要有定位的结果作为依据
***************************************/
int StdVel(ONEpoch* e, SPPResult* r)
{

	//初始化部分局部变量
	vector<double> Lg(e->nGvalid), Mg(e->nGvalid), Ng(e->nGvalid);
	vector<double> Lb(e->nBvalid), Mb(e->nBvalid), Nb(e->nBvalid);
	vector<double> w(e->nGvalid + e->nBvalid);
	for (int i = 1; i < e->nGvalid + 1; i++)//GPS
	{
		double dist = sqrt((r->X - e->gInfo[i].X) * (r->X - e->gInfo[i].X) + (r->Y - e->gInfo[i].Y) * (r->Y - e->gInfo[i].Y) + (r->Z - e->gInfo[i].Z) * (r->Z - e->gInfo[i].Z));
		double dX = r->X - e->gInfo[i].X;
		double dY = r->Y - e->gInfo[i].Y;
		double dZ = r->Z - e->gInfo[i].Z;
		Lg[i - 1] = dX / dist;
		Mg[i - 1] = dY / dist;
		Ng[i - 1] = dZ / dist;
		double wx = -Lg[i - 1] * e->gInfo[i].Vx;
		double wy = -Mg[i - 1] * e->gInfo[i].Vy;
		double wz = -Ng[i - 1] * e->gInfo[i].Vz;
		double rdot = wx + wy + wz;
		double D1 = e->gObs[i].raw[0].dopp * C_Light / FG1_GPS;
		double cdtr = e->gInfo[i].dtr * C_Light;
		w[i - 1] = D1 - rdot + cdtr;
	}
	for (int i = 1; i < e->nBvalid + 1; i++)
	{
		//这是系数的初值
		double dist = sqrt((r->X - e->bInfo[i].X) * (r->X - e->bInfo[i].X) + (r->Y - e->bInfo[i].Y) * (r->Y - e->bInfo[i].Y) + (r->Z - e->bInfo[i].Z) * (r->Z - e->bInfo[i].Z));
		Lb[i - 1] = (r->X - e->bInfo[i].X) / dist;
		Mb[i - 1] = (r->Y - e->bInfo[i].Y) / dist;
		Nb[i - 1] = (r->Z - e->bInfo[i].Z) / dist;
		double wx = -Lb[i - 1] * e->bInfo[i].Vx;
		double wy = -Mb[i - 1] * e->bInfo[i].Vy;
		double wz = -Nb[i - 1] * e->bInfo[i].Vz;
		double D1 = e->bObs[i].raw[0].dopp * C_Light / FB1_BDS;
		double rdot = wx + wy + wz;
		double cdtr = e->bInfo[i].dtr * C_Light;
		w[i + e->nGvalid - 1] = D1 - rdot + cdtr;
	}
	vector<double> Bb(4 * (e->nGvalid + e->nBvalid));

	//分别给BDS和GPS数据赋值系数矩阵
	for (int r = 0; r < e->nGvalid; r++)
	{
		Bb[4 * r] = Lg[r]; Bb[4 * r + 1] = Mg[r]; Bb[4 * r + 2] = Ng[r];
		Bb[4 * r + 3] = 1;
	}
	for (int r = e->nGvalid; r < e->nBvalid + e->nGvalid; r++)
	{
		int s = r - e->nGvalid;
		Bb[4 * r] = Lb[s]; Bb[4 * r + 1] = Mb[s]; Bb[4 * r + 2] = Nb[s];
		Bb[4 * r + 3] = 1;
	}
	CMatrix B(Bb.data(), e->nBvalid + e->nGvalid, 4);
	CMatrix W(w.data(), e->nBvalid + e->nGvalid);
	CMatrix Nbb = B.transpose() * B;
	CMatrix VX = Nbb.inv() * B.transpose() * W;
	CMatrix V = B * VX + (-W);
	double Sig0 = sqrt((V.transpose() * V).value[0] / (e->nBvalid + e->nGvalid - 4));
	CMatrix Q = (B.transpose() * B).inv();
	r->Vx = VX.value[0];
	r->Vy = VX.value[1];
	r->Vz = VX.value[2];
	//for (int r = 0; r < VX.length; r++)
	//	cout << VX.value[r] << endl;
	return 1;
	//	else return 0;
}

void StdResOut(SPPResult* r)
{
	if (r->Sig0 > 100.0)
	{
		cout << "历元解算失败！" << endl; return;
	}
	UTCTime ut = MJD2UTC(GPST2MJD(r->Gtime));
	if (int(ut.Second + 0.1) == 60)
	{
		ut.Minute += 1;
		if (ut.Minute == 60) { ut.Hour += 1; ut.Minute = 0; }
		ut.Second = 0;
	}
	cout << setw(4) << ut.Year << " " << ut.Month << " " << ut.Day << " " << setw(2) << ut.Hour << " " << setw(2) << ut.Minute << " " << setw(2) << int(ut.Second + 0.1) << "  ";
	BLH blh = r->blh;
	cout << setiosflags(ios::fixed) << setprecision(3) << setw(10) << r->X << "  ";
	cout << setw(12) << r->Y << "  ";
	cout << setw(12) << r->Z << "  ";
	cout << setw(12) << setprecision(8) << blh.b * rou << "  ";
	cout << setw(12) << blh.l * rou << "  ";
	cout << setprecision(3) << setw(6) << r->Vx << "  ";
	cout << setw(6) << r->Vy << "  ";
	cout << setw(6) << r->Vz << "  ";
	cout << setprecision(7) << setw(9) << r->Sig0 << "   ";
	cout << setw(6) << r->PDOP << endl;
}

void FileResOut(ofstream& fout, SPPResult* r)
{
	if (r->Sig0 > 100) { return; }
	UTCTime ut = MJD2UTC(GPST2MJD(r->Gtime));
	fout << ut.Year << " " << ut.Month << " " << ut.Day << "," << ut.Hour << " " << ut.Minute << " " << int(ut.Second + 0.1) << ",";
	fout << setiosflags(ios::fixed) << setprecision(3) << r->X << "," << r->Y << "," << r->Z << ",";
	fout << setprecision(7) << r->blh.b * rou << "," << r->blh.l * rou << ",";
	fout << setprecision(3) << r->Vx << "," << r->Vy << "," << r->Vz << ",";
	fout << setprecision(5) << r->Sig0 << "," << r->PDOP << ",";
	fout << r->nGvalid << "," << r->nBvalid << endl;
}