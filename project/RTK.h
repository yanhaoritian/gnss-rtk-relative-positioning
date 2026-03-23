#pragma once
#include"decoding.h"
#include"lambda.h"
#include"sockets.h"
/*RTKObs结构体，
用于直接存储某一系统双频单差观测值*/
struct RTKObs
{
	NAVSYS sys = NUL;				//卫星系统
	GPSInfo GBase, GRover;
	BDSInfo BBase, BRover;			//两个测站GPS和BDS的卫星位置和速度
	unsigned short PRN = 0;			//PRN号
	double Adr1 = 0.0;
	double Psr1 = 0.0;
	double Adr2 = 0.0;
	double Psr2 = 0.0;					//双频伪距、相位站间单差观测值
	double AngH = 0.0;					//卫星高度角
	double sig2 = 0.0;						//通过高度角（或信噪比）为该单差观测赋予的先验单位权方差
	bool Ref = false;						//参考星标记
	double matA[3] = { 0,0,0 };		//站间单差系数矩阵的值
};
/*
RTK历元双差结构体
可储存某一系统中所有的双差观测值
*/
struct RTKEpoch
{
	NAVSYS sys = NUL;						//卫星系统
	unsigned short Refnum = 0;				//参考星PRN号和在所有观测值中的序号
	unsigned short RefPRN = 0;
	double AdrDiff1[MAXBDSPRN];
	double PsrDiff1[MAXBDSPRN];
	double AdrDiff2[MAXBDSPRN];
	double PsrDiff2[MAXBDSPRN];	//双差伪距、相位观测值
	double W[MAXBDSPRN];				//卫地距双差
};


int PPKFile(FILE* FBase, FILE* FRover, unsigned char* buffBase, unsigned char* buffRover, ONEpoch* EBase, ONEpoch* ERover);
int ReadMsgHeader(FILE* F, unsigned char* buff, ONEpoch* e, NavTime* Gtime);
int FindSatsGPS(ONEpoch* A, ONEpoch* B, RTKObs* C);
int FindSatsBDS(ONEpoch* A, ONEpoch* B, RTKObs* C);
int FindRef(RTKObs* C, int num);
int CalDoubleDiff(RTKObs* C, RTKEpoch* E, int num, double WaveLength1, double WaveLength2);
int GetA(XYZ Sat, XYZ SatRef, XYZ Obs, double axyz[3]);
CMatrix GetBMatrix(XYZ x1, XYZ x2, RTKObs* C, RTKEpoch* E, int num);
CMatrix GetQMatrix(RTKObs* C, int num);
int CalEpoch(NAVSYS sys, XYZ xBase, XYZ xRover, RTKObs* C, RTKEpoch* E, int num);
int DectectDoubleTrue(RTKObs* C, int num);
int CheckTime(FILE* FBase, FILE* FRover, unsigned char* buffBase, unsigned char* buffRover, ONEpoch* EBase, ONEpoch* ERover,
	ObsData* gBase1, ObsData* bBase1, ObsData* gRover1, ObsData* bRover1,
	SPPResult* rBase, SPPResult* rRover, NavTime BaseTime, NavTime RoverTime);//检测时间是否对其的函数，如果没有对齐，则继续读数
int RTKSockets(CfgInfo CgBase, CfgInfo CgRover);//网口读取的程序

int RTK2SysLs(XYZ xBase, XYZ* xRover, RTKObs* C1, RTKObs* C2, RTKEpoch* E1, RTKEpoch* E2,
	int num1, int num2, double WaveLength1, double WaveLength2, double WaveLength3, double WaveLength4, ofstream& F, CMatrix* DeltaX, CMatrix* Dxx, double* RMS, double* RDOP, double* Sigma0, double* Ratio);