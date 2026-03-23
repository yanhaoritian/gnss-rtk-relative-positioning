#pragma once
#include<stdio.h>
#include<string.h>
#include<stdint.h>
#include"Matrix.h"
#include"time.h"
#include"coordinate.h"
using namespace std;
//这是使用常量的定义
/* Physical parameters of the Earth, Sun and Moon  */
#define R_WGS84  6378137.0          /* Radius Earth [m]; WGS-84 ,WGS84的椭球长半径 */
#define F_WGS84  (1.0/298.257223563)  /* Flattening; WGS-84  ，WGS-84的椭球扁率 */
#define Omega_WGS 7.2921151467e-5   /*[rad/s], the earth rotation rate ，地球自转角速度*/
#define GM_Earth   398600.5e+9     /* [m^3/s^2]; WGS-84 ，WGS84的地球常数*/
#define GM_JGM3   398600.4415e+9     /* [m^3/s^2]; JGM3  */

//这是CGCS-2000椭球的各项参数，用于北斗系统
/* Physical parameters of the Earth, Sun and Moon  */
#define R_CGS2K  6378137.0          /* Radius Earth [m]; CGCS2000  */
#define F_CGS2K  1.0/298.257222101  /* Flattening; CGCS2000   */
#define Omega_BDS 7.2921150e-5      /*[rad/s], the earth rotation rate */
#define GM_BDS   398600.4418e+9     /* [m^3/s^2]; CGCS2000  */
#define ID_GPS	0							/*GPS的信号类型*/
#define ID_BDS	4							/*BDS的信号类型*/

#define C_Light 299792458		/*定义光速[m/s]*/
/* some constants about GPS satellite signal */
#define  FG1_GPS  1575.42E6             /* L1信号频率 */
#define  FG2_GPS  1227.60E6             /* L2信号频率 */
#define  FG12R    (77/60.0)				/* FG1_Freq/FG2_Freq */
#define  FG12R2   (5929/3600.0)
#define  WL1_GPS  (C_Light/FG1_GPS)	/*L1和L2的波长*/
#define  WL2_GPS  (C_Light/FG2_GPS)

/* some constants about Compass satellite signal */
#define  FB1_BDS  1561.098E6               /* B1信号的基准频率 */
#define  FB2_BDS  1207.140E6               /* B2信号的基准频率 */
#define  FB3_BDS  1268.520E6               /* B3信号的基准频率 */

#define  FC12R    (FG1_CPS/FG2_CPS)       /* FG1_CPS/FG2_CPS */
#define  FC12R2   (FC12R*FC12R)                /* FG1_CPS^2/FG2_CPS^2 */
#define  FC13R    (FG1_CPS/FG3_CPS)       /* FG1_CPS^2/FG3_CPS^2 */
#define  FC13R2   (FC13R*FC13R)
#define  WL1_BDS  (C_Light/FB1_BDS)	/*B1,B2,B3的波长*/
#define  WL2_BDS  (C_Light/FB2_BDS)
#define  WL3_BDS  (C_Light/FB3_BDS)

#define GPST_BDT  14         /* GPS时与北斗时的差值[s] */
#define MAXCHANNUM 36
#define MAXSATNUM  64
#define MAXGPSPRN  32
#define MAXOBSTYPENUM 9
#define MAXGEOPRN  5         /* 最大的GEO卫星号 */ //不完善，BDS-3的GEO被排除在外了
#define MAXBDSPRN 63
#define MAXRAWLEN   40960

//定义解码数据的宏
#define OEM4SYNC1       0xAA    /* oem7/6/4 message start sync code 1 */
#define OEM4SYNC2       0x44    /* oem7/6/4 message start sync code 2 */
#define OEM4SYNC3       0x12    /* oem7/6/4 message start sync code 3 */
#define OEM4HLEN        28      /* oem7/6/4 message header length (bytes) */

#define POLYCRC32 0xEDB88320u	/*CRC32 ploynomial*/
#define ID_BDS_EPH 1696				//北斗星历
#define ID_RANGE 43						//伪距和相位观测值文件
#define ID_GPS_EPH 7						//GPS卫星星历

/*******************
定义两个枚举变量，包括导航系统变量NAVSYS和数据频段FREQ
*********************/
enum NAVSYS { NUL = -1, GPS = 1, BDS = 2 };
enum FREQ { L1 = 1, L2 = 2 };

/*********************************
GPS星历结构体
成员变量：
PRN--卫星的PRN号
tow-星历发布的周内秒
health--卫星健康状况
IODE--卫星的数据龄期
week--GPS周数
toe--星历参考时刻（卫星钟）
A--轨道长半轴，ecc--离心率，dN--平均角速度
M0--toe时刻平近点角
ome--近地点角距
cuc，cus，crc，crs，cic，cis--调谐改正项
toc--卫星钟计算所用的参考时间
I0--轨道倾角的初始值
dI--轨道倾角的变化率
ome0--升交点赤经，dome--升交点赤经变化率
iodc--卫星钟参数发布
tgd--估计的信号硬件延迟
a0,a1,a2--钟差钟漂钟速
N--改正后的平均角速度
URA--用户数据方差
****************************************/
class GPSEph
{
public:
	unsigned long PRN = 0;
	double tow = 0.0;
	unsigned long health = 1;
	unsigned long IODE1 = 0, IODE2 = 0;
	unsigned long week = 0, zweek = 0;
	double toe = 0.0;
	double A = 0.0;
	double dN = 0.0;
	double M0 = 0.0;
	double ecc = 0.0;
	double ome = 0.0;
	double cuc = 0.0; double cus = 0.0;
	double crc = 0.0; double crs = 0.0;
	double cic = 0.0; double cis = 0.0;
	double I0 = 0.0; double dI = 0.0;
	double ome0 = 0; double dome = 0.0;
	unsigned long iodc = 0;
	double toc = 0;
	double tgd = 0;
	double a0 = 0, a1 = 0, a2 = 0;
	int AS = 0;
	double N = 0;
	double URA = 0;
	NavTime Toe, Toc;
};


/*********************************
BDS星历结构体
与GPS不同的参数有：
sqrtA--半长轴的平方根
tgd1，tgd2--硬件延迟的两个参考值
AODC--数据龄期
其余定义相同
****************************************/
class BDSEph
{
public:
	unsigned long PRN = 0;
	unsigned long week = 0;
	double URA = 0;
	unsigned long health = 1;
	double tgd1 = 0; double tgd2 = 0;
	unsigned long AODC = 0;
	double toc = 0.0;
	double a0 = 0, a1 = 0, a2 = 0;
	unsigned long AODE = 0;
	double toe = 0;
	double sqrtA = 0;
	double dN = 0.0;
	double M0 = 0.0;
	double ecc = 0.0;
	double ome = 0.0;
	double cuc = 0.0; double cus = 0.0;
	double crc = 0.0; double crs = 0.0;
	double cic = 0.0; double cis = 0.0;
	double I0 = 0.0; double dI = 0.0;
	double Ome0 = 0; double dOme = 0.0;
	unsigned long iodc = 0;
	NavTime Toe, Toc;
};

/************************************
存储原始观测值的结构体RANGE
存储的数据有：
PRN--卫星的PRN号
psr,psrS--测量得到的伪距和其标准差
adr,adrS--测量的载波相位和标准差，注意载波相位以周为单位
dopp--多普勒效应测量值
ratio,即C/No--载噪比
locktime--信号锁定时间
track--跟踪状态，根据它能读出信号所属的系统Sys和信号类型Type
Sys,Type--从跟踪状态中读取的系统、信号类型
health--该观测数据的健康状况
**************************************/
struct RANGE
{
	unsigned short PRN = 0;
	double psr = 0;
	float psrS = 0;
	double adr = 0;
	float adrS = 0;
	float dopp = 0;
	float ratio = 0;
	float locktime = 0;
	unsigned long track = 0;
	int Sys = 0;
	int Type = 0;
	bool health = false;
};

/****************************************
BDS与GPS卫星信息结构体，用于存储当前某一时刻卫星的位置和速度
PRN--卫星的PRN，Trop--对流层改正
X,Y,Z--卫星在ECEF坐标系的坐标
Vx,Vy,Vz--卫星的速度
dt,dtr--卫星在该时刻的钟差与钟漂
t--当前的参考时刻（只存储周内秒）
ang--卫星的高度角（如果未进行高度角解算这一项不会纳入其中）
BDS还具有额外的Tgd1
****************************************/
struct BDSInfo
{
	int PRN = 0; double Trop = 0;
	double tgd1 = 0;
	double X = 0, Y = 0, Z = 0;
	double Vx = 0, Vy = 0, Vz = 0;
	double dt = 0; double dtr = 0;
	double t = 0;
	double angH = 0;
	bool flag = true;
	//BDSInfo() {
	//	PRN = 0; Trop = 0;
	//	X = 0, Y = 0, Z = 0;
	//	Vx = 0, Vy = 0, Vz = 0;
	//	dt = 0; dtr = 0;
	//	t = 0; angH = 0;
	//	flag = true;
	//}
};
struct GPSInfo
{
	int PRN = 0; double Trop = 0;
	double X = 0, Y = 0, Z = 0;
	double Vx = 0, Vy = 0, Vz = 0;
	double dt = 0; double dtr = 0;
	double t = 0;
	double angH = 0;
	bool flag = true;
};

//一个历元的观测值的结构体
struct ObsData
{
	NAVSYS sys = NUL;
	int Valid = -1;							//观测值的健康状态，默认为-1，刚刚锁定为0，正常锁定为1
	RANGE raw[2];						//存储原始观测数据
	double MW = 0.0;
	double GF = 0.0;
	double MW_Smooth = 0;
	int N_Locked = 0;
	double angH = 0;
};

/********************************************
存储一个历元观测的结构体，都是双系统并行
包含：
卫星星历信息GPS/BDSEph
观测值信息ObsData g/bObs,可用的卫星数目nG/Bvalid
卫星在观测时的位置和速度信息 BDS/GPSInfo
观测时刻GTime（用GPS时）
*********************************************/
struct ONEpoch
{
	GPSEph Gpse[MAXGPSPRN + 1]; bool gEphValid = false;
	BDSEph Bdse[MAXBDSPRN + 1]; bool bEphValid = false;
	ObsData gObs[MAXGPSPRN + 1]; int nGvalid = 0;
	ObsData bObs[MAXBDSPRN + 1]; int nBvalid = 0;
	GPSInfo gInfo[MAXGPSPRN + 1];
	BDSInfo bInfo[MAXBDSPRN + 1];
	NavTime Gtime = NavTime(0, 0);
};

struct SPPResult
{
	double X = 0, Y = 0, Z = 0;
	double PDOP;
	double Sig0;
	double dtG, dtB;
	BLH blh;
	double Vx, Vy, Vz;
	int nGvalid = 0, nBvalid = 0;
	NavTime Gtime;
};
/*******************
下列函数是二进制数据格式转换函数，可以根据数据格式将它读作不同的数据
*********************/
#define U1(p) (*((uint8_t *)(p)))	//unsigned char
#define I1(p) (*((int8_t  *)(p)))		//signed char
static uint16_t U2(uint8_t* p) { uint16_t u; memcpy(&u, p, 2); return u; }	//unsigned short
static uint32_t U4(uint8_t* p) { uint32_t u; memcpy(&u, p, 4); return u; }	//unsigned int
static int32_t  I4(uint8_t* p) { int32_t  i; memcpy(&i, p, 4); return i; }			//int
static float    R4(uint8_t* p) { float    r; memcpy(&r, p, 4); return r; }			//float
static double   R8(uint8_t* p) { double   r; memcpy(&r, p, 8); return r; }		//double

/********************************
解码相关函数的一次性声明
*********************************/
unsigned int check_crc32(unsigned char* buff, int len);									//CRC校验
int DecodeBinMsg(FILE* FObs);																//读取二进制文件数据
int ReadGPSEph(unsigned char* buff, ONEpoch* e);									//读取GPS星历
int ReadBDSEph(unsigned char* buff, ONEpoch* e);									//读取BDS星历
int ReadRange(unsigned char* buff, ObsData* g, ObsData* b);										//读取原始观测值数据
int DetectOutlier(ObsData* g, ObsData* b, ObsData* g1, ObsData* b1, int Nobs);		//对原始观测值数据进行粗差探测并且进行整理
void ObsToEpoch(ObsData* g, ObsData* b, ONEpoch* e);						//将观测历元转换为观测数据
int GPSCal(GPSEph E, GPSInfo* I, double t);												//计算GPS卫星的位置和速度
int BDSCal(BDSEph E, BDSInfo* I, double t);												//计算BDS卫星的位置和速度
void GetSatInfo(ONEpoch* e);
double GIonFree(double a, double b);
double BIonFree(double a, double b);
int CheckEpoch(ONEpoch* e);
int StdPos(ONEpoch* e, SPPResult* sppr);
int StdVel(ONEpoch* e, SPPResult* r);
void StdResOut(SPPResult* r);
void FileResOut(ofstream& fout, SPPResult* r);
int DetectAngH(XYZ* xyz, ONEpoch* e, double S);
double CalDist(GPSInfo e, CMatrix X);
double CalDist(BDSInfo e, CMatrix X);
double CalDist(XYZ p, XYZ q);
double CalNorm3(CMatrix C);
int LeastSquareIteration(ONEpoch* e, CMatrix& B, CMatrix& W, CMatrix& XR0, CMatrix& d_X);
int PsoteriorOutlier(ONEpoch* e, CMatrix V, double Sig0);
//int DecodeOem719Msg(unsigned char buff[], int lenr, int& lenrem, ONEpoch* e, ofstream &fout);
int input_oem6(unsigned char Buff[], int& d, ONEpoch* e, ObsData* g, ObsData* b, ObsData* g1, ObsData* b1, SPPResult* r, ofstream& fout);
int SolveMsg(int msgID, unsigned char* buff, ObsData* g, ObsData* b, ObsData* g1, ObsData* b1, ONEpoch* e, SPPResult* r);