#include"decoding.h"
#include"RTK.h"
#include<cmath>
using namespace std;


/*
	GPS和BDS是两个系统，所以采用两套寻找共视卫星观测值的策略
	找到历元间的共视卫星，并且找出参考卫星的PRN号码
	注意基准站是参数A，流动站是参数B
	此处接受并且计算的信息有：
		基准站和流动站的卫星信息
		卫星的PRN号
		站间单差观测值
	返回值是两个接收机有多少共视卫星
	注意计算方法是流动站减去基准站
*/
int FindSatsGPS(ONEpoch* A, ONEpoch* B, RTKObs* C)
{
	int num = 0;
	for (int i = 1; i <= A->nGvalid; i++)
	{
		if (A->gObs[i].sys != GPS)continue;
		if (A->gObs[i].Valid < 0)continue;
		unsigned short PRN = A->gObs[i].raw[0].PRN;
		int k = 0; bool val = false;
		//if(B->gObs)
		//复制GPS观测值,计算单差观测值
		for (k = 1; k <= B->nGvalid; k++)
		{
			if (B->gObs[k].raw[0].PRN == PRN)
			{
				val = true;
				break;
			}
		}
		if (val == true)
		{
			if (B->gObs[k].angH < 15.0 / rou)
				continue;
			C[num].sys = GPS;
			C[num].GBase = A->gInfo[i];//C[num].BBase = A->bInfo[i];
			C[num].GRover = B->gInfo[k];//C[num].BRover = B->bInfo[i];
			C[num].PRN = A->gObs[i].raw[0].PRN;
			C[num].Adr1 = -A->gObs[i].raw[0].adr + B->gObs[k].raw[0].adr;
			C[num].Adr2 = -A->gObs[i].raw[1].adr + B->gObs[k].raw[1].adr;
			C[num].Psr1 = -A->gObs[i].raw[0].psr + B->gObs[k].raw[0].psr;
			C[num].Psr2 = -A->gObs[i].raw[1].psr + B->gObs[k].raw[1].psr;
			C[num].AngH = B->gObs[k].angH;
			//if (C[num].AngH < 10.0/rou)
			//	continue;
			C[num].sig2 = (16 + 9.0 / (sin(C[num].AngH) * sin(C[num].AngH))) * 1e-6;
			//C[num].sig2 = 9e-6;
			num++;
		}
	}
	return num;
}

int FindSatsBDS(ONEpoch* A, ONEpoch* B, RTKObs* C)
{
	int num = 0;
	for (int i = 1; i <= A->nBvalid; i++)
	{
		if (A->bObs[i].sys != BDS)continue;
		if (A->bObs[i].Valid < 0)continue;
		//复制BDS观测值,计算单差
		unsigned short PRN = A->bObs[i].raw[0].PRN;
		int k = 0; bool val = false;
		for (k = 1; k <= B->nBvalid; k++)
		{
			if (B->bObs[k].raw[0].PRN == PRN)
			{
				val = true;
				break;
			}
		}
		if (val == true)
		{
			if (B->bObs[k].angH < 15.0 / rou)
				continue;
			C[num].sys = BDS;
			C[num].PRN = A->bObs[i].raw[0].PRN;
			//C[num].GBase = A->gInfo[i]; 
			C[num].BBase = A->bInfo[i]; C[num].BRover = B->bInfo[k];
			C[num].Adr1 = -A->bObs[i].raw[0].adr + B->bObs[k].raw[0].adr;
			C[num].Adr2 = -A->bObs[i].raw[1].adr + B->bObs[k].raw[1].adr;
			C[num].Psr1 = -A->bObs[i].raw[0].psr + B->bObs[k].raw[0].psr;
			C[num].Psr2 = -A->bObs[i].raw[1].psr + B->bObs[k].raw[1].psr;
			C[num].AngH = B->bObs[k].angH;

			C[num].sig2 = (16 + 9.0 / (sin(C[num].AngH) * sin(C[num].AngH))) * 1e-6;
			//C[num].sig2 = 9e-6;
			num++;
		}
	}
	return num;
}
/*
寻找参考星的函数FindRef
	参数：
	RTK观测值数组C
	可用共视卫星数量num（由之前寻找共视卫星的）
*/
int FindRef(RTKObs* C, int num)
{
	if (num <= 1)return 0;
	for (int i = 0; i < num; i++)C[i].Ref = false;
	int PRef = 0; double AngTemp = 0.0;
	int ntemp = 0;
	for (int i = 0; i < num; i++)
	{
		if (AngTemp < C[i].AngH)//&&C[i].PRN!=39)
		{
			AngTemp = C[i].AngH;
			PRef = C[i].PRN;
			ntemp = i;
		}
		//C[i].Ref = true;
	}
	C[ntemp].Ref = true;
	return PRef;
}

/*
计算双差观测值
输入参数：
	存有站间单差观测值的数组C
	RTK单历元双差结构体E
	共视卫星数目num
功能：对RTK双差结构体进行赋值更新
注意双差是观测卫星减去参考卫星，即参考性、基准站都是位于减号后
*/
int CalDoubleDiff(RTKObs* C, RTKEpoch* E, int num, double WaveLength1, double WaveLength2)
{
	E->sys = C[0].sys;
	int Refn = 0;
	for (; Refn < num; Refn++)//寻找参考星在结构体中的位置
	{
		if (C[Refn].Ref == true)
		{
			E->Refnum = Refn;
			E->RefPRN = C[Refn].PRN;
			break;
		}
	}
	if (Refn < 0 || Refn>num)return 0;
	RTKObs Temp = C[Refn]; int j = 0;
	//C[Refn]
	E->sys = C[j].sys;
	for (int i = 0; i < num; i++)
	{
		if (i == Refn)continue;
		E->AdrDiff1[j] = (C[i].Adr1 - C[Refn].Adr1) * WaveLength1;
		E->PsrDiff1[j] = C[i].Psr1 - C[Refn].Psr1;
		E->AdrDiff2[j] = (C[i].Adr2 - C[Refn].Adr2) * WaveLength2;
		E->PsrDiff2[j] = C[i].Psr2 - C[Refn].Psr2;
		j++;

	}
	return 1;
}
/*
完成双差解算之后，再计算常数项（距离之差）为最小二乘做好准备
参数：
	sys-系统
	xBase,xRover - 坐标位置参数（基准站和流动站）
	C - RTK观测值数组
	E - RTK 单历元结构体
	num - 实际测得的有效共视卫星
*/
int CalEpoch(NAVSYS sys, XYZ xBase, XYZ xRover, RTKObs* C, RTKEpoch* E, int num)
{
	int r = E->Refnum; int j = 0;
	XYZ xRefB(0, 0, 0), xRefR(0, 0, 0); double distBR = 0.0, distRR = 0.0, distBS = 0.0, distRS = 0.0;
	if (sys == NUL)return -1;
	if (sys == GPS)
	{
		xRefB = XYZ(C[r].GBase.X, C[r].GBase.Y, C[r].GBase.Z);
		xRefR = XYZ(C[r].GRover.X, C[r].GRover.Y, C[r].GRover.Z);
		distBR = CalDist(xRefB, xBase);
		distRR = CalDist(xRefR, xRover);
		for (int i = 0; i < num; i++)
		{
			if (i == r)continue;
			XYZ  xSatB(C[i].GBase.X, C[i].GBase.Y, C[i].GBase.Z);
			XYZ  xSatR(C[i].GRover.X, C[i].GRover.Y, C[i].GRover.Z);
			distRS = CalDist(xSatR, xRover);
			distBS = CalDist(xSatB, xBase);
			double Ll = distRS - distBS - distRR + distBR;//常数项，几何距离之差
			E->W[j] = Ll; E->W[j + num - 1] = Ll;//双频观测有时需要重复两次
			j++;
		}
	}
	if (sys == BDS)
	{
		xRefB = XYZ(C[r].BBase.X, C[r].BBase.Y, C[r].BBase.Z);
		xRefR = XYZ(C[r].BRover.X, C[r].BRover.Y, C[r].BRover.Z);
		distBR = CalDist(xRefB, xBase);
		distRR = CalDist(xRefR, xRover);
		for (int i = 0; i < num; i++)
		{
			if (i == r)continue;
			XYZ  xSatB(C[i].BBase.X, C[i].BBase.Y, C[i].BBase.Z);
			XYZ  xSatR(C[i].BRover.X, C[i].BRover.Y, C[i].BRover.Z);
			distRS = CalDist(xSatR, xRover);
			distBS = CalDist(xSatB, xBase);
			double Ll = distRS - distBS - distRR + distBR;//常数项，几何距离之差
			E->W[j] = Ll; E->W[j + num - 1] = Ll;
			j++;
		}
	}
	return 1;
}

int DectectDoubleTrue(RTKObs* C, int num)
{
	int r = 0;
	for (int i = 0; i < num; i++)
	{
		if (C[i].Ref == true)r++;
	}
	return r;
}

