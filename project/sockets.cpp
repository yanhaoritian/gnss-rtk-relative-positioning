#include"sockets.h"
#include"decoding.h"
#include<iostream>
#include<fstream>
using namespace std;


bool OpenSocket(SOCKET& sock, const char IP[], const unsigned short Port)
{
	WSADATA wsaData;
	SOCKADDR_IN addrSrv = {};

	if (!WSAStartup(MAKEWORD(1, 1), &wsaData))
	{
		if ((sock = socket(AF_INET, SOCK_STREAM, 0)) != INVALID_SOCKET)
		{
			addrSrv.sin_addr.S_un.S_addr = inet_addr(IP);
			addrSrv.sin_family = AF_INET;
			addrSrv.sin_port = htons(Port);
			if (connect(sock, (SOCKADDR*)&addrSrv, sizeof(SOCKADDR)) == 0)
			{
				return true;
			}
			closesocket(sock);
		}
	}
	return false;
}

void CloseSocket(SOCKET& sock)
{
	closesocket(sock);
	WSACleanup();
}

/* input oem6 raw data from iostream ------------------------------------------
input oem6����������Զ�ȡʵʱ��OEM719��Ϣ
���������
ʵʱ������������Buff
��һ�λ�����δ��������ݳ���d
����Ԫ�۲�ṹ��ONEpoch *e
�۲�ṹ�壺ObsData��GPS��*g����BDS��*b��
��һ��Ԫ�Ĺ۲�ṹ�壺*g1��*b1
������ļ�������ofstream &fout
*-----------------------------------------------------------------------------*/
int input_oem6(unsigned char Buff[], int& d, ONEpoch* e, ObsData* g, ObsData* b, ObsData* g1, ObsData* b1, SPPResult* r, ofstream& fout)
{
	int i, j, len, Nrange, val;
	double tow;
	int msgType, week, msgID;
	unsigned char TempBuff[MAXRAWLEN];
	int saved = d;
	i = 0;
	val = 0;

	int count = 0;

	while (1)
	{
		for (; i < d - 2; i++) //ͬ��
		{
			if (Buff[i] == OEM4SYNC1 && Buff[i + 1] == OEM4SYNC2 && Buff[i + 2] == OEM4SYNC3)break;
			//�����������Ϣͷ����ô�ͽ���������
		}
		if (i + OEM4HLEN >= d)		break;										//���ʣ������ݲ����洢һ����Ϣͷ

		for (j = 0; j < OEM4HLEN; j++)TempBuff[j] = Buff[i + j];//�洢��Ϣͷ
		len = U2(TempBuff + 8) + OEM4HLEN;

		if ((len + 4 + i) > d || len > MAXRAWLEN)break;			//��Ϣ������������ 


		for (j = OEM4HLEN; j < len + 4; j++)								 //�����Ϣ���壬Ҳ��һ��Tempbuff��ֻ��洢һ����Ϣͷ����Ϣ���壨��CRCУ��λ��
			TempBuff[j] = Buff[i + j];
		msgID = U2(TempBuff + 4);


		/* check crc32 */
		if (check_crc32(TempBuff, len) != U4(TempBuff + len))
		{
			i += len + 4;
			continue;
		}
		msgType = (U1(TempBuff + 6) >> 4) & 0x03;
		week = U2(TempBuff + 14);
		tow = U4(TempBuff + 16) * 0.001;

		if (msgType != 0)
			continue; /* message type: 0=binary,1=ascii */
		//cout << "ID" << msgID << endl;
		switch (msgID)											//������Ϣ����ѡ����
		{
		case ID_BDS_EPH:
			ReadBDSEph(TempBuff, e);
			break;
		case ID_GPS_EPH:
			ReadGPSEph(TempBuff, e);
			break;
		case ID_RANGE:
			e->Gtime = NavTime(week, tow);
			val = 1;	//α��ȹ۲��� Prn�����ֽ�
			Nrange = ReadRange(TempBuff, g, b);
			DetectOutlier(g, b, g1, b1, Nrange);
			ObsToEpoch(g, b, e);
			if (CheckEpoch(e))//�����������ͨ��
			{
				StdPos(e, r);
				StdVel(e, r);
				//StdResOut(&r);
				FileResOut(fout, r);
				//cout << e->Gtime.WeekSec << "  ";
				//cout << MJD2UTC(GPST2MJD(e->Gtime)).Second << "  " << endl;
				for (int i = 1; i < MAXGPSPRN + 1; i++) { g1[i] = g[i]; e->gInfo[i].PRN = 0; }
				for (int i = 1; i < MAXBDSPRN + 1; i++) { b1[i] = b[i]; e->bInfo[i].PRN = 0; }
			}


			break;
		default:
			break;
		}
		i += len + 4;

		if (val == 1)  //����ɹ�
			break;
	}

	for (j = 0; j < saved - i; j++)
		Buff[j] = Buff[i + j];

	d = j; //����󣬻�����ʣ�����δ������ֽ���
	if (d < 0)d = 0;
	return val;
}