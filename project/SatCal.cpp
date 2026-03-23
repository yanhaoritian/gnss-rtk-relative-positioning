#include"decoding.h"
#include<cmath>
#include"Matrix.h"
using namespace std;
/********************************
GPS卫星位置速度计算函数，将原始星历数据转换成含有卫星速度和位置信息的结构体
主要是计算卫星位置和速度、改正钟差和钟速
每次只确定单卫星的位置和速度
**********************************/
int GPSCal(GPSEph E, GPSInfo* G, double t)
{

	//含有摄动项的二体问题解决卫星的位置
	double n0 = sqrt(GM_Earth / pow(E.A, 3));			//计算平均转速
	double tk = t - E.toe;										    //计算时间差
	//double tkc = t - E.toc;

	if (tk > 302400)	tk -= 604800;
	else if (tk < -302400)	tk += 604800;						//对周数进行必要的判断

	double n = n0 + E.dN;											//改正平均角速度
	double Mk = E.M0 + n * tk;
	double Ek0 = 0, Ek = Mk;
	double e = E.ecc;
	int j = 0; double dEk = 1;
	while (dEk > pow(10, -12))
	{
		Ek0 = Mk + e * sin(Ek);
		dEk = fabs(Ek - Ek0);
		Ek = Ek0;
		j++;
		if (j > 200)
		{
			cout << "警告：开普勒方程陷入死循环" << endl;
			break;
		}
	}

	double vk = atan2(sqrt(1 - e * e) * sin(Ek) / (1 - e * cos(Ek)), (cos(Ek) - e) / (1 - e * cos(Ek)));
	double Pk = vk + E.ome;
	double duk = E.cus * sin(2 * Pk) + E.cuc * cos(2 * Pk);
	double drk = E.crs * sin(2 * Pk) + E.crc * cos(2 * Pk);
	double dik = E.cis * sin(2 * Pk) + E.cic * cos(2 * Pk);
	double uk = Pk + duk;
	double rk = E.A * (1 - e * cos(Ek)) + drk;
	double ik = E.I0 + dik + E.dI * tk;
	double xk = rk * cos(uk);
	double yk = rk * sin(uk);
	double Omek = E.ome0 + (E.dome - Omega_WGS) * tk - Omega_WGS * E.toe;
	G->X = xk * cos(Omek) - yk * cos(ik) * sin(Omek);
	G->Y = xk * sin(Omek) + yk * cos(ik) * cos(Omek);
	G->Z = yk * sin(ik);

	//卫星的钟差改正
	double tRelate = -4.442807633e-10 * e * sqrt(E.A) * sin(Ek);//相对论改正
	double dtsv = E.a0 + E.a1 * (t - E.toc) + E.a2 * (t - E.toc) * (t - E.toc) + tRelate;// -E.tgd;
	G->dt = dtsv;

	//计算卫星速度
	double dEkd = n / (1 - e * cos(Ek));
	double dPhk = sqrt((1 + e) / (1 - e)) * (cos(vk / 2) / cos(Ek / 2)) * (cos(vk / 2) / cos(Ek / 2)) * dEkd;
	double dukd = 2 * (E.cus * cos(2 * Pk) - E.cuc * sin(2 * Pk)) * dPhk + dPhk;
	double drkd = E.A * e * sin(Ek) * dEkd + 2 * (E.crs * cos(2 * Pk) - E.crc * sin(2 * Pk)) * dPhk;
	double dIkd = E.dI + 2 * (E.cis * cos(2 * Pk) - E.cic * sin(2 * Pk)) * dPhk;
	double dOme = E.dome - Omega_WGS;
	double Rvalue[] = {
		cos(Omek), -sin(Omek) * cos(ik), -(xk * sin(Omek) + yk * cos(Omek) * cos(ik)), yk * sin(Omek) * sin(ik),
		sin(Omek),cos(Omek) * cos(ik), (xk * cos(Omek) - yk * sin(Omek) * cos(ik)), -yk * cos(Omek) * sin(ik),
		0,sin(ik),0,yk * cos(ik)
	};
	CMatrix Rdot(Rvalue, 3, 4);
	double dxk[] = {
		drkd * cos(uk) - rk * dukd * sin(uk),
		drkd * sin(uk) + rk * dukd * cos(uk),
		dOme,
		dIkd };
	CMatrix S(dxk, 4, 1);
	CMatrix X(3, 1);
	X = Rdot * S;
	G->Vx = X.value[0];
	G->Vy = X.value[1];
	G->Vz = X.value[2];

	double dtR = -4.442807633e-10 * e * sqrt(E.A) * Ek * dEkd;
	double ddtsv = E.a1 + 2 * E.a2 * (t - E.toc) + dtR;
	G->dtr = ddtsv;
	return 1;
}

/****************************************
//北斗卫星位置和速度的计算程序
//BDS里面PRN1-5，59-61都是GEO
北斗要格外注意GEO的解算问题，GPS和BDS的位置算法可以在ICD文档中找到
而北斗的GEO卫星要通过轨道旋转5度在惯性系中解算
****************************************/
int BDSCal(BDSEph E, BDSInfo* B, double t)
{

	double A = E.sqrtA * E.sqrtA;
	//含有摄动项的二体问题解决卫星的位置
	double n0 = sqrt(GM_BDS / pow(A, 3));
	double tk = t - E.toe;
	if (tk > 302400)	tk -= 604800;
	else if (tk < -302400)	tk += 604800;
	double n = n0 + E.dN;
	double Mk = E.M0 + n * tk;
	double Ek0 = 0, Ek = Mk;
	double e = E.ecc;
	int j = 0; double dEk = 1;
	while (dEk > pow(10, -12))
	{
		Ek0 = Mk + e * sin(Ek);
		dEk = abs(Ek - Ek0);
		Ek = Ek0;
		j++;
		if (j > 200)
		{
			cout << "警告：开普勒方程陷入死循环" << endl;
			break;
		}
	}

	double vk = atan2(sqrt(1 - e * e) * sin(Ek) / (1 - e * cos(Ek)), (cos(Ek) - e) / (1 - e * cos(Ek)));
	double Pk = vk + E.ome;
	double duk = E.cus * sin(2 * Pk) + E.cuc * cos(2 * Pk);
	double drk = E.crs * sin(2 * Pk) + E.crc * cos(2 * Pk);
	double dik = E.cis * sin(2 * Pk) + E.cic * cos(2 * Pk);
	double uk = Pk + duk;
	double rk = A * (1 - e * cos(Ek)) + drk;
	double ik = E.I0 + dik + E.dI * tk;
	double xk = rk * cos(uk);
	double yk = rk * sin(uk);
	if ((E.PRN > 5 && E.PRN < 59) || E.PRN > 61)				//不是GEO的普通卫星（MEO和IGSO）
	{
		double Omek = E.Ome0 + (E.dOme - Omega_BDS) * tk - Omega_BDS * E.toe;
		B->X = xk * cos(Omek) - yk * cos(ik) * sin(Omek);
		B->Y = xk * sin(Omek) + yk * cos(ik) * cos(Omek);
		B->Z = yk * sin(ik);

		//卫星的钟差改正
		double tRelate = -4.442807633e-10 * e * sqrt(A) * sin(Ek);//相对论改正
		double dtsv = E.a0 + E.a1 * (t - E.toc) + E.a2 * (t - E.toc) * (t - E.toc) + tRelate;// -E.tgd1;
		B->dt = dtsv;

		//计算卫星速度
		double dEkd = n / (1 - e * cos(Ek));
		double dPhk = sqrt((1 + e) / (1 - e)) * (cos(vk / 2) / cos(Ek / 2)) * (cos(vk / 2) / cos(Ek / 2)) * dEkd;
		double dukd = 2 * (E.cus * cos(2 * Pk) - E.cuc * sin(2 * Pk)) * dPhk + dPhk;
		double drkd = A * e * sin(Ek) * dEkd + 2 * (E.crs * cos(2 * Pk) - E.crc * sin(2 * Pk)) * dPhk;
		double dIkd = E.dI + 2 * (E.cis * cos(2 * Pk) - E.cic * sin(2 * Pk)) * dPhk;
		double dOme = E.dOme - Omega_BDS;
		double Rvalue[] = {
			cos(Omek), -sin(Omek) * cos(ik), -(xk * sin(Omek) + yk * cos(Omek) * cos(ik)), yk * sin(Omek) * sin(ik),
			sin(Omek),cos(Omek) * cos(ik), (xk * cos(Omek) - yk * sin(Omek) * cos(ik)), yk * cos(Omek) * sin(ik),
			0,sin(ik),0,yk * cos(ik)
		};
		CMatrix Rdot(Rvalue, 3, 4);
		double dxk[] = {
			drkd * cos(uk) - rk * dukd * sin(uk),
			drkd * sin(uk) + rk * dukd * cos(uk),
			dOme,
			dIkd };
		CMatrix S(dxk, 4, 1);
		CMatrix X = Rdot * S;
		B->Vx = X.value[0];
		B->Vy = X.value[1];
		B->Vz = X.value[2];

		double dtR = -4.442807633e-10 * e * sqrt(A) * Ek * dEkd;
		double ddtsv = E.a1 + 2 * E.a2 * (t - E.toc) + dtR;
		B->dtr = ddtsv;

	}
	else//如果是GEO
	{
		double Omek = E.Ome0 + E.dOme * tk - Omega_BDS * E.toe;//计算历元升交点经度（惯性系）
		double Xgk = xk * cos(Omek) - yk * cos(ik) * sin(Omek);
		double Ygk = xk * sin(Omek) + yk * cos(ik) * cos(Omek);
		double Zgk = yk * sin(ik);
		double xyzgk[] = { Xgk,Ygk,Zgk };
		CMatrix gk(xyzgk, 3);
		double sinp = sin(Omega_BDS * tk), cosp = cos(Omega_BDS * tk);
		double rz[] =
		{ cosp,sinp,0,
			-sinp,cosp,0,
			0,0,1 };
		double sin5 = sin(-5.0 / rou);
		double cos5 = cos(-5.0 / rou);
		double rx[] =
		{ 1,0,0,
			0,cos5,sin5,
			0,-sin5,cos5 };//旋转-5度
		CMatrix Rz(rz, 3, 3); CMatrix Rx(rx, 3, 3);
		CMatrix XK = Rz * Rx * gk;
		B->X = XK.value[0];
		B->Y = XK.value[1];
		B->Z = XK.value[2];

		//卫星的钟差改正
		double tRelate = -4.442807633e-10 * e * sqrt(A) * sin(Ek);//相对论改正
		double dtsv = E.a0 + E.a1 * (t - E.toc) + E.a2 * (t - E.toc) * (t - E.toc) + tRelate;// -E.tgd1;
		B->dt = dtsv;
		//卫星的速度修正
		double dEkd = n / (1 - e * cos(Ek));
		double dPhk = sqrt((1 + e) / (1 - e)) * (cos(vk / 2) / cos(Ek / 2)) * (cos(vk / 2) / cos(Ek / 2)) * dEkd;
		double dukd = 2 * (E.cus * cos(2 * Pk) - E.cuc * sin(2 * Pk)) * dPhk + dPhk;
		double drkd = A * e * sin(Ek) * dEkd + 2 * (E.crs * cos(2 * Pk) - E.crc * sin(2 * Pk)) * dPhk;
		double dIkd = E.dI + 2 * (E.cis * cos(2 * Pk) - E.cic * sin(2 * Pk)) * dPhk;
		double dOme = E.dOme;
		double Rvalue[] = {
			cos(Omek), -sin(Omek) * cos(ik), -(xk * sin(Omek) + yk * cos(Omek) * cos(ik)), yk * sin(Omek) * sin(ik),
			sin(Omek),cos(Omek) * cos(ik), (xk * cos(Omek) - yk * sin(Omek) * cos(ik)), yk * cos(Omek) * sin(ik),
			0,sin(ik),0,yk * cos(ik)
		};
		CMatrix Rdot(Rvalue, 3, 4);
		double dxk[] = {
			drkd * cos(uk) - rk * dukd * sin(uk),
			drkd * sin(uk) + rk * dukd * cos(uk),
			dOme,
			dIkd };
		CMatrix S(dxk, 4, 1);
		CMatrix X = Rdot * S;

		double thetaE = Omega_BDS * tk;
		double st = sin(thetaE); double ct = cos(thetaE);
		double r1[] = {
			ct,st * cos5,st * sin5,
			-st,ct * cos5,ct * sin5,
			0,-sin5,cos5
		};
		CMatrix R1(r1, 3, 3);
		X = R1 * X;
		double r2[] = {
		Omega_BDS * (-st) * Xgk + Omega_BDS * ct * cos5 * Ygk + Omega_BDS * ct * sin5 * Zgk,
		Omega_BDS * (-ct) * Xgk + Omega_BDS * (-st) * cos5 * Ygk + Omega_BDS * (-st) * sin5 * Zgk,
		0 };
		CMatrix R2(r2, 3);
		CMatrix VX = R2 + X;
		double dtR = -4.442807633e-10 * e * sqrt(A) * Ek * dEkd;
		double ddtsv = E.a1 + 2 * E.a2 * (t - E.toc) + dtR;
		B->dtr = ddtsv;
		B->Vx = VX.value[0];
		B->Vy = VX.value[1];
		B->Vz = VX.value[2];

	}

	return 1;
}

/***********************************************
一个观测历元数据中自行解算卫星位置和速度的参数
基于上述两个函数进行解算
最终将位置、速度、钟差存入ONEpoch结构体的info
同时把信号发射时刻对应的周内秒也存储进去
此处暂时不考虑迭代,仅使用一次的值
*************************************************/
void GetSatInfo(ONEpoch* e)
{
	for (int i = 1; i < e->nGvalid + 1; i++)
	{
		//计算信号发射时刻的位置和速度
		int PRN = e->gObs[i].raw[0].PRN;
		//if (e->Gpse[PRN].PRN == 0)continue;
		e->gInfo[i].PRN = PRN;
		double t = e->Gtime.WeekSec - e->gObs[i].raw[0].psr / C_Light;
		GPSCal(e->Gpse[PRN], &e->gInfo[i], t);
		t = t - e->gInfo[i].dt;
		GPSCal(e->Gpse[PRN], &e->gInfo[i], t);
		e->gInfo[i].t = t;
		double Dt = e->gObs[i].raw[0].psr / C_Light + e->gInfo[i].dt;//信号传播的时间
		//进行地球自转改正
		double wt = Omega_WGS * Dt;
		double rr[] =
		{
			cos(wt),sin(wt),0,
			-sin(wt),cos(wt),0
			,0,0,1
		};
		CMatrix Rr(rr, 3, 3);
		//Rr = Eye(3);
		double xs0[] = { e->gInfo[i].X,e->gInfo[i].Y,e->gInfo[i].Z };
		CMatrix Xs0(xs0, 3);
		CMatrix Xs = Rr * Xs0;
		e->gInfo[i].X = Xs.value[0];
		e->gInfo[i].Y = Xs.value[1];
		e->gInfo[i].Z = Xs.value[2];
	}
	for (int i = 1; i < e->nBvalid + 1; i++)
	{
		//计算信号发射时刻的卫星位置和速度

		int PRN = e->bObs[i].raw[0].PRN;
		e->bInfo[i].PRN = PRN;
		//if (e->Bdse[PRN].PRN == 0)continue;
		e->bInfo[i].tgd1 = e->Bdse[PRN].tgd1;
		NavTime Btime = GPS2BDS(e->Gtime);
		double t_pse = e->bObs[i].raw[0].psr / C_Light;
		double t = Btime.WeekSec - t_pse;
		BDSCal(e->Bdse[PRN], &e->bInfo[i], t);
		t = t - e->bInfo[i].dt;
		BDSCal(e->Bdse[PRN], &e->bInfo[i], t);
		e->bInfo[i].t = t;
		double Dt = e->bObs[i].raw[0].psr / C_Light + e->bInfo[i].dt;
		//进行地球自转改正
		double wt = Omega_BDS * Dt;
		double rr[] = {
			cos(wt),sin(wt),0,
			-sin(wt),cos(wt),0,
			0,0,1 };
		CMatrix Rr(rr, 3, 3);
		double xs0[] = { e->bInfo[i].X,e->bInfo[i].Y,e->bInfo[i].Z };
		CMatrix Xs0(xs0, 3);
		CMatrix Xs = Rr * Xs0;
		e->bInfo[i].X = Xs.value[0];
		e->bInfo[i].Y = Xs.value[1];
		e->bInfo[i].Z = Xs.value[2];

	}

}