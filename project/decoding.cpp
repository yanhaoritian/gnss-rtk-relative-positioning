#include<windows.h>
#include"sockets.h"
#include"decoding.h"
#include<iostream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<stdlib.h>

using namespace std;

//CRC-32位校验码
unsigned int check_crc32(unsigned char* buff, int len)
{
	int i, j;
	unsigned int crc = 0;
	for (i = 0; i < len; i++)
	{
		crc ^= buff[i];
		for (j = 0; j < 8; j++)
		{
			if (crc & 1)crc = (crc >> 1) ^ POLYCRC32;
			else crc >>= 1;
		}
	}
	return crc;
}

int ReadMsgHeader(FILE* F, unsigned char* buff, ONEpoch* e, NavTime* Gtime)
{
	if (fread(buff + 3, 1, 25, F) < 1)return -2;
	int num = 4;
	unsigned short msgID = U2(buff + num);									num += 2;
	char msgType = I1(buff + num);													num += 1;
	int type = (msgType >> 4) & (0x3);
	char portAdd = U1(buff + num);													num += 1;
	unsigned short msgLength = U2(buff + num);							num += 2;
	unsigned short Sequence = U2(buff + num);								num += 2;
	unsigned char IdleTime = U1(buff + num);									num += 1; num += 1;
	unsigned short week = U2(buff + num);										num += 2;
	unsigned long ms = U4(buff + num);											num += 4;
	NavTime T = NavTime(week, ms / 1000.0);
	e->Gtime = T;
	*Gtime = NavTime(T.Week, T.WeekSec);
	num += 8;
	if (fread(buff + num, 1, msgLength + 4, F) < 1)return -2;
	unsigned int CRC0 = U4(buff + num + msgLength);
	unsigned int CRC = check_crc32(buff, msgLength + 28);//CRC
	if (CRC0 != CRC)
	{
		cerr << "CRC Check Failed!" << endl;
		return -1;
	}
	return msgID;

}
/******************************************
读取二进制消息的函数，这是封装外围的函数之一，只需要传递文件指针参数即可
它阅读每一条二进制消息，并将二进制消息转换为结构体
********************************************/
int DecodeBinMsg(FILE* FObs)
{
	unsigned char buff[MAXRAWLEN] = { 0 };					//创建等同于最大原始数据数的空数组
	//GPSEph* gpse = new GPSEph[MAXGPSPRN+1];
	//BDSEph* bdse = new BDSEph[MAXBDSPRN+1];	
	SPPResult r;
	ONEpoch e;
	ObsData g1[MAXGPSPRN + 1];							//存储上一个历元的观测值，创建大小等同于卫星的最大PRN的数组
	ObsData b1[MAXBDSPRN + 1];

	ofstream fout("StandardPointPositioningRes.csv");
	fout << "日期,时间,ECEF-X/m,ECEF-Y/m,ECEF-Z/m,ECEF-B/m,ECEF-L/m,Vx(m/s),Vy(m/s),Vz(m/s),Sigma0/m,PDOP/m,GPS卫星数目,BDS卫星数目" << endl;
	//初始化部分必需数据
	while (1)
	{
		if (fread(buff + 2, 1, 1, FObs) < 1)return 1;	//读不到数据的时候就返回
		if (buff[0] == OEM4SYNC1 && buff[1] == OEM4SYNC2 && buff[2] == OEM4SYNC3)
		{

			if (fread(buff + 3, 1, 25, FObs) < 1)return 1;
			int num = 4;
			unsigned short msgID = U2(buff + num);									num += 2;
			char msgType = I1(buff + num);													num += 1;
			int type = (msgType >> 4) & (0x3);
			char portAdd = U1(buff + num);													num += 1;
			unsigned short msgLength = U2(buff + num);							num += 2;
			unsigned short Sequence = U2(buff + num);								num += 2;
			unsigned char IdleTime = U1(buff + num);									num += 1; num += 1;
			unsigned short week = U2(buff + num);										num += 2;
			unsigned long ms = U4(buff + num);											num += 4;
			NavTime Gt = NavTime(week, ms / 1000.0);
			e.Gtime = Gt;
			num += 8;
			if (fread(buff + num, 1, msgLength + 4, FObs) < 1)return -1;
			unsigned int CRC0 = U4(buff + num + msgLength);
			unsigned int CRC = check_crc32(buff, msgLength + 28);//CRC
			if (CRC0 != CRC)
			{
				cerr << "CRC Check Failed!" << endl;
				continue;
			}
			//RANGE *range = new RANGE[MAXSATNUM * 10];
			int Nrange = 0;
			ObsData g[MAXGPSPRN + 1];								//创建大小等同于卫星的最大PRN的数组
			ObsData b[MAXBDSPRN + 1];


			switch (msgID)
			{
			case ID_BDS_EPH:
				ReadBDSEph(buff, &e);
				break;
			case ID_GPS_EPH:
				ReadGPSEph(buff, &e);
				break;
			case ID_RANGE:
				Nrange = ReadRange(buff, g, b);
				DetectOutlier(g, b, g1, b1, Nrange);
				ObsToEpoch(g, b, &e);
				if (CheckEpoch(&e))
				{
					StdPos(&e, &r);
					StdVel(&e, &r);
					StdResOut(&r);
					FileResOut(fout, &r);
				}
				for (int i = 1; i < MAXGPSPRN + 1; i++) { g1[i] = g[i]; e.gInfo[i].PRN = 0; }
				for (int i = 1; i < MAXBDSPRN + 1; i++) { b1[i] = b[i]; e.bInfo[i].PRN = 0; }
				break;
			default: break;
			}
		}
		else//如果找不到匹配的字符，则将buff移位
		{
			buff[0] = buff[1];
			buff[1] = buff[2];
		}
	}
	fout.close();
	return 0;
}

int ReadGPSEph(unsigned char* buff, ONEpoch* e)
{
	if (buff == NULL)return -1;
	buff += 28;
	int nsat = U4(buff);						buff += 4;
	if (nsat <= 0 || nsat > MAXGPSPRN) return -1;
	GPSEph gpse;
	gpse.PRN = nsat;
	gpse.tow = R8(buff);			buff += 8;
	gpse.health = U4(buff);		buff += 4;
	gpse.IODE1 = U4(buff);		buff += 4;
	gpse.IODE2 = U4(buff);		buff += 4;
	gpse.week = U4(buff);			buff += 4;
	gpse.zweek = U4(buff);		buff += 4;
	gpse.toe = R8(buff);			buff += 8;
	gpse.A = R8(buff);				buff += 8;
	gpse.dN = R8(buff);			buff += 8;
	gpse.M0 = R8(buff);			buff += 8;
	gpse.ecc = R8(buff);			buff += 8;
	gpse.ome = R8(buff);			buff += 8;
	gpse.cuc = R8(buff);			buff += 8;
	gpse.cus = R8(buff);			buff += 8;
	gpse.crc = R8(buff);				buff += 8;
	gpse.crs = R8(buff);				buff += 8;
	gpse.cic = R8(buff);				buff += 8;
	gpse.cis = R8(buff);				buff += 8;
	gpse.I0 = R8(buff);				buff += 8;
	gpse.dI = R8(buff);				buff += 8;
	gpse.ome0 = R8(buff);			buff += 8;
	gpse.dome = R8(buff);			buff += 8;
	gpse.iodc = U4(buff);			buff += 4;
	gpse.toc = R8(buff);			buff += 8;
	gpse.tgd = R8(buff);			buff += 8;
	gpse.a0 = R8(buff);				buff += 8;
	gpse.a1 = R8(buff);				buff += 8;
	gpse.a2 = R8(buff);				buff += 8;
	gpse.AS = U4(buff);			buff += 4;
	gpse.N = R8(buff);				buff += 8;
	gpse.URA = R8(buff);			buff += 8;
	unsigned long CRC32 = U4(buff);	buff += 4;
	gpse.Toc = NavTime(static_cast<unsigned short>(gpse.week), gpse.toc);
	gpse.Toe = NavTime(static_cast<unsigned short>(gpse.week), gpse.toe);
	e->Gpse[nsat] = gpse;
	e->gEphValid = true;
	//if (gpse[nsat].health != 0)
	//{
	//	e->Gpse[nsat].PRN = 0;
	//	e->Gpse[nsat].health = 0;
	//}
	return 0;
}

int ReadBDSEph(unsigned char* buff, ONEpoch* e)
{
	if (buff == NULL)return -1;
	buff += 28;
	unsigned long nsat = U4(buff);
	if (nsat <= 0 || nsat > MAXBDSPRN) return -1;
	BDSEph bdse;
	bdse.PRN = nsat;				buff += 4;
	bdse.week = U4(buff);			buff += 4;
	bdse.URA = R8(buff);			buff += 8;
	bdse.health = U4(buff);		buff += 4;
	bdse.tgd1 = R8(buff);			buff += 8;
	bdse.tgd2 = R8(buff);			buff += 8;
	bdse.AODC = U4(buff);		buff += 4;
	bdse.toc = U4(buff);			buff += 4;
	bdse.a0 = R8(buff);				buff += 8;
	bdse.a1 = R8(buff);				buff += 8;
	bdse.a2 = R8(buff);				buff += 8;
	bdse.AODE = U4(buff);		buff += 4;
	bdse.toe = U4(buff);			buff += 4;
	bdse.sqrtA = R8(buff);			buff += 8;
	bdse.ecc = R8(buff);			buff += 8;
	bdse.ome = R8(buff);			buff += 8;
	bdse.dN = R8(buff);			buff += 8;
	bdse.M0 = R8(buff);			buff += 8;
	bdse.Ome0 = R8(buff);		buff += 8;
	bdse.dOme = R8(buff);		buff += 8;
	bdse.I0 = R8(buff);				buff += 8;
	bdse.dI = R8(buff);				buff += 8;
	bdse.cuc = R8(buff);			buff += 8;
	bdse.cus = R8(buff);			buff += 8;
	bdse.crc = R8(buff);			buff += 8;
	bdse.crs = R8(buff);				buff += 8;
	bdse.cic = R8(buff);				buff += 8;
	bdse.cis = R8(buff);				buff += 8;
	unsigned long CRC32 = U4(buff);	buff += 4;
	bdse.Toc = NavTime(static_cast<unsigned short>(bdse.week), bdse.toc);
	bdse.Toe = NavTime(static_cast<unsigned short>(bdse.week), bdse.toe);
	e->Bdse[nsat] = bdse;
	if (bdse.health == 1)
	{
		e->Bdse[nsat].PRN = 0;
		e->Bdse[nsat].health = 1;
	}
	e->bEphValid = true;
	return 0;
}

int ReadRange(unsigned char* buff, ObsData* g, ObsData* b)
{
	if (buff == NULL)return -1;
	buff += 28;
	unsigned long Nobs = U4(buff);				buff += 4;
	if (Nobs > MAXSATNUM * 10) Nobs = MAXSATNUM * 10;
	RANGE r[MAXSATNUM * 10] = {};
	for (unsigned int i = 0; i < Nobs; i++)
	{
		unsigned short PRN = U2(buff);
		//if (PRN == 33) { buff += 44; continue; }
		r[i].PRN = PRN;				buff += 2; buff += 2;
		r[i].psr = R8(buff);			buff += 8;
		r[i].psrS = R4(buff);			buff += 4;
		r[i].adr = -R8(buff);			buff += 8;
		r[i].adrS = R4(buff);			buff += 4;
		r[i].dopp = -R4(buff);		buff += 4;
		r[i].ratio = R4(buff);			buff += 4;
		r[i].locktime = R4(buff);	buff += 4;
		r[i].track = U4(buff);		buff += 4;
		r[i].Sys = (r[i].track >> 0x10) & (0x7);
		r[i].Type = (r[i].track >> 0x15) & (0x1F);
		int phlock = (r[i].track >> (0xA)) & (0x1);			//探测载波相位和伪距的失锁情况
		int rlock = (r[i].track >> (0xC)) & (0x1);
		if (r[i].ratio <= 0 || r[i].locktime <= 0)r[i].health = false;
		else if (rlock == 0)r[i].health = false;
		else r[i].health = true;
	}
	for (unsigned long i = 0; i < Nobs; i++)
	{
		int Freq = -1;
		unsigned short PRN = r[i].PRN;
		if (r[i].Sys == ID_GPS)//GPS
		{
			switch (r[i].Type)
			{
			case 0:  Freq = 0; break;   // L1C/A
			case 9:  Freq = 1; break;   // L2P(Y),semi-codeless
			default: Freq = -1; break;
			}
			if (Freq < 0 || Freq>1)continue;
			if (Freq == 0 && r[i].ratio < 35)continue;
			if (Freq == 1 && r[i].ratio < 20)continue;
			g[PRN].sys = GPS;
			g[PRN].raw[Freq] = r[i];
		}
		else if (r[i].Sys == ID_BDS)//BDS
		{
			switch (r[i].Type)
			{
			case 0: Freq = 0; break;   // B1I D1
			case 2: Freq = 1; break;   // B3I D1
			case 4: Freq = 0; break;   // B1I D2
			case 6: Freq = 1; break;   // B3I D2
			default: Freq = -1; break;
			}
			if (Freq < 0 || Freq>1)continue;
			if (Freq == 0 && r[i].ratio < 35)continue;
			//if (Freq == 1 && r[i].ratio < 25)continue;
			b[PRN].sys = BDS;
			b[PRN].raw[Freq] = r[i];
		}
		else continue;
	}
	//	g[2]; b[3];
	return Nobs;
}

int SolveMsg(int msgID, unsigned char* buff, ObsData* g, ObsData* b, ObsData* g1, ObsData* b1, ONEpoch* e, SPPResult* r)
{
	int Nrange = 0;
	switch (msgID)
	{
	case ID_GPS_EPH:
		ReadGPSEph(buff, e);
		break;
	case ID_BDS_EPH:
		ReadBDSEph(buff, e);
		break;
	case ID_RANGE:
		Nrange = ReadRange(buff, g, b);
		DetectOutlier(g, b, g1, b1, Nrange);
		ObsToEpoch(g, b, e);
		if (CheckEpoch(e))
		{
			StdPos(e, r);
			StdVel(e, r);
		}
		else return -1;
		for (int i = 1; i < MAXGPSPRN + 1; i++) g1[i] = g[i];
		for (int i = 1; i < MAXBDSPRN + 1; i++) b1[i] = b[i];
		break;
	default:
		return -1;
	}
	return 1;
}
