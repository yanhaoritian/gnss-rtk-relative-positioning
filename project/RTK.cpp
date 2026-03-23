#include"decoding.h"
#include"RTK.h"
#include<iostream>
#include <fstream>
#include<iomanip>

using namespace std;

int CheckTime(FILE* FBase, FILE* FRover, unsigned char* buffBase, unsigned char* buffRover, ONEpoch* EBase, ONEpoch* ERover,
	ObsData* gBase1, ObsData* bBase1, ObsData* gRover1, ObsData* bRover1,
	SPPResult* rBase, SPPResult* rRover, NavTime BaseTime, NavTime RoverTime)
{
	while (1)
	{
		if (BaseTime - RoverTime > 0.01)//基准站的时间更靠后，此时应该继续读流动站的数据以进行时间对齐
		{
			if (fread(buffRover + 2, 1, 1, FRover) < 1)return -1;	//读不到数据的时候就返回
			if (buffRover[0] == OEM4SYNC1 && buffRover[1] == OEM4SYNC2 && buffRover[2] == OEM4SYNC3)
			{
				int Nrange = 0;
				ObsData gRover[MAXGPSPRN + 1]; ObsData bRover[MAXBDSPRN + 1]; 								//创建大小等同于卫星的最大PRN的数组
				int msgID = ReadMsgHeader(FRover, buffRover, ERover, &RoverTime);
				SolveMsg(msgID, buffRover, gRover, bRover, gRover1, bRover1, ERover, rRover);
				if (msgID == ID_RANGE)break;
			}
			else//如果找不到消息头，就令buffer移位
			{
				buffRover[0] = buffRover[1];
				buffRover[1] = buffRover[2];
			}
		}
		else if (BaseTime - RoverTime < -0.01)//流动站的时间更靠后，此时应该读取基准站的数据进行时间对齐
		{
			if (fread(buffBase + 2, 1, 1, FBase) < 1)return -1;	//读不到数据的时候就返回
			if (buffBase[0] == OEM4SYNC1 && buffBase[1] == OEM4SYNC2 && buffBase[2] == OEM4SYNC3)
			{
				ObsData gBase[MAXGPSPRN + 1]; ObsData bBase[MAXBDSPRN + 1]; 								//创建大小等同于卫星的最大PRN的数组
				int msgID = ReadMsgHeader(FBase, buffBase, EBase, &BaseTime);
				SolveMsg(msgID, buffBase, gBase, bBase, gBase1, bBase1, EBase, rBase);
				if (msgID == ID_RANGE)break;
			}
			else//如果找不到消息头，就令buffer移位
			{
				buffBase[0] = buffBase[1];
				buffBase[1] = buffBase[2];
			}
		}
		else return 1;
	}

	return 1;
}

int PPKFile(FILE* FBase, FILE* FRover, unsigned char* buffBase, unsigned char* buffRover, ONEpoch* EBase, ONEpoch* ERover)
{
	ObsData gBase1[MAXGPSPRN + 1]; ObsData bBase1[MAXBDSPRN + 1]; 								//创建大小等同于卫星的最大PRN的数组
	ObsData gRover1[MAXGPSPRN + 1]; ObsData bRover1[MAXBDSPRN + 1];
	SPPResult rBase, rRover; RTKObs RgObs[MAXGPSPRN + 1], RbObs[MAXBDSPRN + 1];
	ofstream f("RTKRes.csv"); NavTime BaseTime, RoverTime;
	f << "GPS周," << "GPS周秒," << "X-ECEF," << "Y-ECEF," << "Z-ECEF," << "N," << "E," << "U," << "Sigx," << "Sigy," << "Sigz," << "RDOP," << "PDOP," << "RMS," << "AR_STATE," << "RATIO" << endl;
	ofstream F("RTKrecord.txt");
	while (1)
	{
		int SPPFlagBase = -1, SPPFlagRover = -1;
		CheckTime(FBase, FRover, buffBase, buffRover, EBase, ERover, gBase1, bBase1, gRover1, bRover1, &rBase, &rRover, BaseTime, RoverTime);
		while (1)
		{
			if (fread(buffBase + 2, 1, 1, FBase) < 1)return -1;	//读不到数据的时候就返回
			if (buffBase[0] == OEM4SYNC1 && buffBase[1] == OEM4SYNC2 && buffBase[2] == OEM4SYNC3)
			{
				ObsData gBase[MAXGPSPRN + 1]; ObsData bBase[MAXBDSPRN + 1]; 								//创建大小等同于卫星的最大PRN的数组
				int msgID = ReadMsgHeader(FBase, buffBase, EBase, &BaseTime);
				SPPFlagBase = SolveMsg(msgID, buffBase, gBase, bBase, gBase1, bBase1, EBase, &rBase);
				if (msgID == ID_RANGE)break;
			}
			else//如果找不到消息头，就令buffer移位
			{
				buffBase[0] = buffBase[1];
				buffBase[1] = buffBase[2];
			}
		}
		while (1)
		{
			if (fread(buffRover + 2, 1, 1, FRover) < 1)return -1;	//读不到数据的时候就返回
			if (buffRover[0] == OEM4SYNC1 && buffRover[1] == OEM4SYNC2 && buffRover[2] == OEM4SYNC3)
			{
				int Nrange = 0;
				ObsData gRover[MAXGPSPRN + 1]; ObsData bRover[MAXBDSPRN + 1]; 								//创建大小等同于卫星的最大PRN的数组
				int msgID = ReadMsgHeader(FRover, buffRover, ERover, &RoverTime);
				SPPFlagRover = SolveMsg(msgID, buffRover, gRover, bRover, gRover1, bRover1, ERover, &rRover);
				if (msgID == ID_RANGE)break;
			}
			else//如果找不到消息头，就令buffer移位
			{
				buffRover[0] = buffRover[1];
				buffRover[1] = buffRover[2];
			}
		}

		if (SPPFlagBase > 0 && SPPFlagRover > 0)
		{
			NavTime GTime = ERover->Gtime;
			int nGPS = FindSatsGPS(EBase, ERover, RgObs);
			int nBDS = FindSatsBDS(EBase, ERover, RbObs);
			int gRef = FindRef(RgObs, nGPS);
			int bRef = FindRef(RbObs, nBDS);
			for (int i = 1; i < MAXGPSPRN + 1; i++) { EBase->gInfo[i].PRN = 0; ERover->gInfo[i].PRN = 0; }
			for (int i = 1; i < MAXBDSPRN + 1; i++) { EBase->bInfo[i].PRN = 0; ERover->bInfo[i].PRN = 0; }
			for (int i = 0; i < nGPS; i++)
				cout << "G" << RgObs[i].PRN << " ";
			//cout << endl;
			for (int i = 0; i < nBDS; i++)
				cout << "C" << RbObs[i].PRN << " ";
			cout << endl;
			cout << nGPS << " " << gRef << " " << nBDS << " " << bRef << endl;
			XYZ xyzBase(-2267804.5263, 5009342.3723, 3220991.8632);
			XYZ xyzRover(rRover.X, rRover.Y, rRover.Z);
			//cout << rRover.X + 2267804.5263 << "," << rRover.Y - 5009342.3723<<"," << rRover.Z - 3220991.8632 << endl;
			RTKEpoch Eb, Eg;
			if (nBDS > 1 && nGPS > 1)
			{
				//for(int i=0;i<MAXBDSPRN;)
				CalDoubleDiff(RbObs, &Eb, nBDS, WL1_BDS, WL3_BDS);
				//CalDoubleDiff(RbObs, &Eb, nBDS, WL1_BDS);
				CalEpoch(BDS, xyzBase, xyzRover, RbObs, &Eb, nBDS);
				CalDoubleDiff(RgObs, &Eg, nGPS, WL1_GPS, WL2_GPS);
				CalEpoch(GPS, xyzBase, xyzRover, RgObs, &Eg, nGPS);
				if (DectectDoubleTrue(RbObs, nBDS) > 1)
					cout << "B" << endl;
				if (DectectDoubleTrue(RgObs, nGPS) > 1)
					cout << "G" << endl;
				CMatrix dX, Dxx; double Ratio = 0; double RMS = 1000, Sig0 = 1000, RDOP = 1000;
				
				F << setw(4) << GTime.Week << "\t" << GTime.WeekSec << endl;
				//dX = RTKLS(xyzBase, &xyzRover, RgObs, &Eg, nGPS, WL1_GPS);
				int State = RTK2SysLs(xyzBase, &xyzRover, RgObs, RbObs, &Eg, &Eb, nGPS, nBDS, WL1_GPS, WL2_GPS, WL1_BDS, WL3_BDS, F, &dX, &Dxx, &RMS, &RDOP, &Sig0, &Ratio);
				//int State = RTK2SysLs(xyzBase, &xyzRover, RgObs, RbObs, &Eg, &Eb, nGPS, nBDS, WL1_GPS, WL2_GPS, WL1_BDS, F, &dX, &Dxx, &RMS, &RDOP, &Sig0, &Ratio);
				//MatrixShow(dX);
				cout << Ratio << endl;
				//XYZ dx(dX.value[0],dX.value[1], dX.value[2]);
				XYZ Xfixed = xyzRover - xyzBase;
				NEU neu;
				XYZ2NEU(WGS84, &xyzRover, &xyzBase, &neu);
				f << GTime.Week << "," << GTime.WeekSec << ",";
				cout << setw(4) << GTime.Week << "," << GTime.WeekSec << endl;
				cout << setiosflags(ios::fixed) << setprecision(5) << Xfixed.x << "\t" << Xfixed.y << "\t" << Xfixed.z << endl;
				f << setiosflags(ios::fixed) << setprecision(5) << Xfixed.x << "," << Xfixed.y << "," << Xfixed.z << ",";
				f << setiosflags(ios::fixed) << setprecision(5) << neu.n << "," << neu.e << "," << neu.u << ",";
				f << setiosflags(ios::fixed) << setprecision(5) << sqrt(Dxx.find(1, 1)) << "," << sqrt(Dxx.find(2, 2)) << "," << sqrt(Dxx.find(3, 3))
					<< "," << RDOP << "," << sqrt(tr(Dxx) / Sig0) << ",";
				f << RMS << "," << State << "," << Ratio << endl;
				F << setiosflags(ios::fixed) << setprecision(5) << Xfixed.x << "," << Xfixed.y << "," << Xfixed.z << ",";
				F << Ratio << endl;
				F << endl;
			}
		}
	}

	return 1;
}

