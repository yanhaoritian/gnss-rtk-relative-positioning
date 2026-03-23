#include"decoding.h"
#include"RTK.h"
#include"sockets.h"
#include<iostream>
#include<fstream>
#include<iomanip>
#include<cstring>
using namespace std;
namespace {
	constexpr double kTimeAlignThresholdSec = 0.01;
	constexpr int kRecvSleepMs = 900;
	enum class SyncState {
		SyncBoth = 0,
		CatchRover = 1,
		CatchBase = -1
	};

	void AppendRecvBuffer(unsigned char* dst, int& savedLen, const unsigned char* src, int readLen)
	{
		if (readLen <= 0) return;
		if (savedLen < 0 || savedLen + readLen > 2 * MAXRAWLEN) savedLen = 0;
		memcpy(dst + savedLen, src, readLen);
		savedLen += readLen;
	}

	void ResetEpochSatInfo(ONEpoch* base, ONEpoch* rover)
	{
		for (int i = 1; i <= MAXGPSPRN; ++i) { base->gInfo[i].PRN = 0; rover->gInfo[i].PRN = 0; }
		for (int i = 1; i <= MAXBDSPRN; ++i) { base->bInfo[i].PRN = 0; rover->bInfo[i].PRN = 0; }
	}

	void OutputRTKResultLine(ofstream& fres, const NavTime& time, const XYZ& xyzRover, const XYZ& xyzBase, const CMatrix& Dxx, double RDOP, double Sig0, double RMS, int state, double ratio)
	{
		XYZ fixed = xyzRover - xyzBase;
		NEU neu;
		XYZ roverCopy = xyzRover;
		XYZ baseCopy = xyzBase;
		XYZ2NEU(WGS84, &roverCopy, &baseCopy, &neu);
		fres << time.Week << "," << time.WeekSec << ",";
		fres << setiosflags(ios::fixed) << setprecision(5) << fixed.x << "," << fixed.y << "," << fixed.z << ",";
		fres << setiosflags(ios::fixed) << setprecision(5) << neu.n << "," << neu.e << "," << neu.u << ",";
		fres << setiosflags(ios::fixed) << setprecision(5) << sqrt(Dxx.find(1, 1)) << "," << sqrt(Dxx.find(2, 2)) << "," << sqrt(Dxx.find(3, 3))
			<< "," << RDOP << "," << sqrt(tr(Dxx) / Sig0) << ",";
		fres << RMS << "," << state << "," << ratio << endl;
	}
}

int RTKSockets(CfgInfo CgBase, CfgInfo CgRover)
{
	FILE* BaseRecv = fopen("BaseRecv.log", "wb");
	FILE* RoverRecv = fopen("RoveRecv.log", "wb");
	ofstream RoverSppfout("RoverSPPResult.csv");
	RoverSppfout << "日期,时间,ECEF-X/m,ECEF-Y/m,ECEF-Z/m,ECEF-B/m,ECEF-L/m,Vx(m/s),Vy(m/s),Vz(m/s),Sigma0/m,PDOP/m,GPS卫星数目,BDS卫星数目" << endl;
	ofstream BaseSppfout("BaseSPPResult.csv");
	BaseSppfout << "日期,时间,ECEF-X/m,ECEF-Y/m,ECEF-Z/m,ECEF-B/m,ECEF-L/m,Vx(m/s),Vy(m/s),Vz(m/s),Sigma0/m,PDOP/m,GPS卫星数目,BDS卫星数目" << endl;
	ofstream fres("RTKRes.csv");
	ofstream fout("RTKFloat.txt");
	fres << "GPS周," << "GPS周秒," << "X-ECEF," << "Y-ECEF," << "Z-ECEF," << "N," << "E," << "U," << "Sigx," << "Sigy," << "Sigz," << "RDOP," << "PDOP," << "RMS," << "AR_STATE," << "RATIO" << endl;
	WSACleanup();
	//fwrite(buff, sizeof(unsigned char), lenR, Fobs);
	if (OpenSocket(CgBase.Net, CgBase.NetIP, CgBase.NetPort) == false) printf("The ip %s was not opened\n", CgBase.NetIP);
	//ConfigSocket(CgBase.Net);
	if (OpenSocket(CgRover.Net, CgRover.NetIP, CgRover.NetPort) == false) printf("The ip %s was not opened\n", CgRover.NetIP);
	ONEpoch EBase, ERover;
	ObsData gBase1[MAXGPSPRN + 1], bBase1[MAXBDSPRN + 1], gRover1[MAXGPSPRN + 1], bRover1[MAXBDSPRN + 1];
	int lenReadBase = 0, lenSavedBase = 0; int lenReadRover = 0, lenSavedRover = 0;
	unsigned char buffBase[MAXRAWLEN], buffBaseRes[2 * MAXRAWLEN];
	unsigned char buffRover[MAXRAWLEN], buffRoverRes[2 * MAXRAWLEN];
	SyncState state = SyncState::SyncBoth;
	while (1)
	{
		Sleep(kRecvSleepMs);
		int BaseReadFlag, RoverReadFlag;
		SPPResult rBase, rRover;
		//Sleep(300);
		if (state == SyncState::SyncBoth)
		{
			lenReadBase = recv(CgBase.Net, (char*)buffBase, MAXRAWLEN, 0);
			fwrite(buffBase, sizeof(unsigned char), lenReadBase, BaseRecv);
			lenReadRover = recv(CgRover.Net, (char*)buffRover, MAXRAWLEN, 0);
			if (lenReadRover > 1e7 || lenReadRover < 0)lenReadRover = 0;
			fwrite(buffRover, sizeof(unsigned char), lenReadRover, RoverRecv);
		}
		else if (state == SyncState::CatchRover)
		{
			lenReadRover = recv(CgRover.Net, (char*)buffRover, MAXRAWLEN, 0);
			//lenReadBase = 0;
			fwrite(buffRover, sizeof(unsigned char), lenReadRover, RoverRecv);
		}
		else if (state == SyncState::CatchBase)
		{
			lenReadBase = recv(CgBase.Net, (char*)buffBase, MAXRAWLEN, 0);
			//lenReadRover = 0;
			fwrite(buffBase, sizeof(unsigned char), lenReadBase, BaseRecv);
		}
		if (lenReadBase > 0)
			AppendRecvBuffer(buffBaseRes, lenSavedBase, buffBase, lenReadBase);
		if (lenReadRover > 0)
			AppendRecvBuffer(buffRoverRes, lenSavedRover, buffRover, lenReadRover);
		if (lenReadRover > 0 || lenReadBase > 0)
		{
			ObsData gBase[MAXGPSPRN + 1], bBase[MAXBDSPRN + 1];
			ObsData gRover[MAXGPSPRN + 1], bRover[MAXBDSPRN + 1];
			while (1)
			{
				if (state == SyncState::SyncBoth)
				{
					BaseReadFlag = input_oem6(buffBaseRes, lenSavedBase, &EBase, gBase, bBase, gBase1, bBase1, &rBase, BaseSppfout);
					RoverReadFlag = input_oem6(buffRoverRes, lenSavedRover, &ERover, gRover, bRover, gRover1, bRover1, &rRover, RoverSppfout);
					if (lenSavedBase < 0)lenSavedBase = 0;
					if (lenSavedRover < 0)lenSavedRover = 0;
					double DT = EBase.Gtime - ERover.Gtime;
					if (DT > kTimeAlignThresholdSec) { state = SyncState::CatchRover; }// 基准站时间更靠后（Tbase > Trover）
					else if (DT < -kTimeAlignThresholdSec) { state = SyncState::CatchBase; }// 流动站时间更靠后（Tbase < Trover）
					else
					{
						RTKObs RgObs[MAXGPSPRN], RbObs[MAXBDSPRN];
						int nGPS = FindSatsGPS(&EBase, &ERover, RgObs);
						int nBDS = FindSatsBDS(&EBase, &ERover, RbObs);
						int gRef = FindRef(RgObs, nGPS);
						int bRef = FindRef(RbObs, nBDS);
						//cout << EBase.Gtime.Week << "," << EBase.Gtime.WeekSec << endl;
						ResetEpochSatInfo(&EBase, &ERover);
						for (int i = 0; i < nGPS; i++)
							cout << "G" << RgObs[i].PRN << " ";
						for (int i = 0; i < nBDS; i++)
							cout << "C" << RbObs[i].PRN << " ";
						cout << endl;
						cout << nGPS << " " << gRef << " " << nBDS << " " << bRef << endl;
						XYZ xyzBase(-2267804.5263, 5009342.3723, 3220991.8632);
						XYZ xyzRover(rRover.X, rRover.Y, rRover.Z);
						RTKEpoch Eb, Eg;
						if (nBDS > 1 && nGPS > 1)
						{
							CalDoubleDiff(RbObs, &Eb, nBDS, WL1_BDS, WL3_BDS);
							CalEpoch(BDS, xyzBase, xyzRover, RbObs, &Eb, nBDS);
							CalDoubleDiff(RgObs, &Eg, nGPS, WL1_GPS, WL2_GPS);
							CalEpoch(GPS, xyzBase, xyzRover, RgObs, &Eg, nGPS);
							if (DectectDoubleTrue(RbObs, nBDS) > 1)
								cout << "BERR" << endl;
							if (DectectDoubleTrue(RgObs, nGPS) > 1)
								cout << "GERR" << endl;
							CMatrix dX, Dxx; double Ratio = -1; double RMS = 9999, Sig0 = 9999, RDOP = 9999;
							int State = RTK2SysLs(xyzBase, &xyzRover, RgObs, RbObs, &Eg, &Eb, nGPS, nBDS, WL1_GPS, WL2_GPS, WL1_BDS, WL3_BDS, fout, &dX, &Dxx, &RMS, &RDOP, &Sig0, &Ratio);
							XYZ dx = xyzRover - xyzBase;
							cout << dx.x << " " << dx.y << " " << dx.z << endl;
							cout << Ratio << endl;
							OutputRTKResultLine(fres, EBase.Gtime, xyzRover, xyzBase, Dxx, RDOP, Sig0, RMS, State, Ratio);
							cout << endl;
							state = SyncState::SyncBoth;
						}
						break;
					}
				}
				else if (state == SyncState::CatchRover)
				{
					RoverReadFlag = input_oem6(buffRoverRes, lenSavedRover, &ERover, gRover, bRover, gRover1, bRover1, &rRover, fout);
					if (lenSavedRover < 0)lenSavedRover = 0;
					//if (RoverReadFlag == 0)break;
					double DT = EBase.Gtime - ERover.Gtime;
					if (DT > kTimeAlignThresholdSec) { state = SyncState::CatchRover; }// 基准站时间更靠后（Tbase > Trover）
					else if (DT < -kTimeAlignThresholdSec) { state = SyncState::CatchBase; }
					else { state = SyncState::SyncBoth; }
					break;
				}
				else if (state == SyncState::CatchBase)
				{
					BaseReadFlag = input_oem6(buffBaseRes, lenSavedBase, &EBase, gBase, bBase, gBase1, bBase1, &rBase, fout);
					for (int i = 1; i < MAXBDSPRN + 1; i++)ERover.Bdse[i] = EBase.Bdse[i];
					for (int i = 1; i < MAXGPSPRN + 1; i++)ERover.Gpse[i] = EBase.Gpse[i];
					ERover.bEphValid = EBase.bEphValid; ERover.gEphValid = EBase.gEphValid;
					if (lenSavedBase < 0)lenSavedBase = 0;
					//if (BaseReadFlag == 0)break;
					double DT = EBase.Gtime - ERover.Gtime;
					if (DT > kTimeAlignThresholdSec) { state = SyncState::CatchRover; }// 基准站时间更靠后（Tbase > Trover）
					else if (DT < -kTimeAlignThresholdSec) { state = SyncState::CatchBase; }
					else { state = SyncState::SyncBoth; }
					break;
				}
			}
		}
	}
	fclose(BaseRecv); fclose(RoverRecv);
	return 1;
}