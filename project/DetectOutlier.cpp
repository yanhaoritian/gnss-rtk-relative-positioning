#include"decoding.h"
/*MW组合和GF的组合探测*/
double MWgroup(ObsData* Obs, NAVSYS a)
{
	double FL1, FL2, LL1, LL2, MW;
	if (a == GPS)
	{
		FL1 = FG1_GPS;
		FL2 = FG2_GPS;
		LL1 = WL1_GPS;
		LL2 = WL2_GPS;
	}
	else if (a == BDS)
	{
		FL1 = FB1_BDS;
		FL2 = FB3_BDS;
		LL1 = WL1_BDS;
		LL2 = WL3_BDS;
	}
	else return 0;
	MW = (FL1 - FL2) / (FL1 + FL2) * ((Obs->raw[0].psr / LL1) + (Obs->raw[1].psr / LL2)) - (Obs->raw[0].adr - Obs->raw[1].adr);
	Obs->MW = MW;
	return MW;
}
double GFgroup(ObsData* Obs)
{
	double GF = Obs->raw[0].psr - Obs->raw[1].psr;
	Obs->GF = GF;
	return GF;
}
/*******************************************
粗差探测和观测历元结构体的整理函数，如果粗差探测过关，再将其整理进ONEpoch结构体中
******************************************/
int DetectOutlier(ObsData* g, ObsData* b, ObsData* g1, ObsData* b1, int Nobs)
{
	//进行MW和GF的组合观测检验
	for (int i = 1; i < MAXGPSPRN + 1; i++)
	{
		if (g[i].raw[0].health == true && g[i].raw[1].health == true)
		{
			double GF = GFgroup(&g[i]);
			double MW = MWgroup(&g[i], GPS);
			double dGF = 100.0, dMW = 100.0;
			if (g1[i].Valid >= 0)
			{
				g[i].Valid = 1;
				g[i].N_Locked = g1[i].N_Locked + 1;
				dGF = GF - g1[i].GF;
				dMW = MW - g1[i].MW_Smooth;
				g[i].MW_Smooth = (g1[i].MW_Smooth * g1[i].N_Locked + g[i].MW) / g[i].N_Locked;
				if (dGF > 0.05 || dMW > 10)
					g[i].Valid = -1;
			}
			else//刚刚观测到这颗卫星，让它也参与解算
			{
				g[i].N_Locked = 0;
				g[i].Valid = 0;
				g[i].GF = GF;
				g[i].MW = MW;
				g[i].MW_Smooth = MW;
			}
		}
		else
		{
			g[i].Valid = -1;
			g[i].N_Locked = 0;

		}
		//double MW

	}
	for (int i = 1; i < MAXBDSPRN + 1; i++)
	{
		if (b[i].raw[0].health == true && b[i].raw[1].health == true)
		{
			double GF = GFgroup(&b[i]);
			double MW = MWgroup(&b[i], BDS);
			double dGF = 100, dMW = 100;
			if (b1[i].Valid >= 0)
			{
				b[i].Valid = 1;
				b[i].N_Locked = b1[i].N_Locked + 1;
				dGF = GF - b1[i].GF;
				//double 
				dMW = MW - b1[i].MW_Smooth;
				b[i].MW_Smooth = (b1[i].MW_Smooth * b1[i].N_Locked + b[i].MW) / b[i].N_Locked;
				if (dGF > 0.05 || dMW > 10)
					b[i].Valid = -1;
			}
			else//刚刚观测到这颗卫星，让它也参与解算
			{
				b[i].N_Locked = 0;
				b[i].Valid = 0;
				b[i].GF = GF;
				b[i].MW = MW;
				b[i].MW_Smooth = MW;
			}
		}
		else
		{
			b[i].Valid = -1;
			b[i].N_Locked = 0;
		}
	}

	return 1;
}
void ObsToEpoch(ObsData* g, ObsData* b, ONEpoch* e)
{

	int tempg = 1, tempb = 1;
	for (int n = 1; n <= MAXGPSPRN; n++)
	{
		double tk = e->Gtime.WeekSec - e->Gpse[n].toe;
		if (tk > 302400)tk -= 604800;
		else if (tk < -302400)tk += 604800;
		int weekErr = e->Gtime.Week - e->Gpse[n].week;
		if (weekErr > 1 || weekErr < -1)
			continue;
		if (fabs(tk) >= 7500)continue;
		//如果卫星星历和观测值都是健康的话，将其计为可用历元
		if (g[n].Valid >= 0 && e->Gpse[n].PRN != 0)
		{
			e->gObs[tempg] = g[n];
			tempg++;
		}
	}
	e->nGvalid = tempg - 1;
	for (int n = 1; n <= MAXBDSPRN; n++)
	{
		NavTime Bdt = GPS2BDS(e->Gtime);
		int weekErr = Bdt.Week - e->Bdse[n].week;
		if (weekErr > 1 || weekErr < -1)
			continue;
		double tk = Bdt.WeekSec - e->Bdse[n].toe;
		if (tk > 302400)tk -= 604800;
		else if (tk < -302400)tk += 604800;
		if (fabs(tk) >= 3900)continue;
		//如果观测值和星历都健康的话，将其计为可用的历元
		if (b[n].Valid >= 0 && e->Bdse[n].PRN != 0)
		{
			//b[n].Valid = true;
			e->bObs[tempb] = b[n];
			tempb++;
		}

	}
	e->nBvalid = tempb - 1;
}
/**********************************
//计算卫星的高度角进行排除的函数
它用在SPP完成一次解算迭代之后的高度角粗差排除
这里直接对历元卫星信息结构体以及B，W矩阵进行更新
**********************************/
int DetectAngH(XYZ* xyz, ONEpoch* e, double S)
{
	//CMatrix D = B;
	//double X = xyz->x; double y = xyz->y; double z = xyz->z;
	int temp = 0; double ang;
	int nGtemp = e->nGvalid; int nBtemp = e->nBvalid;
	for (int i = 1; i < nGtemp + 1; i++)
	{
		unsigned short PRN = e->gInfo[i].PRN;
		XYZ xyz1(e->gInfo[i].X, e->gInfo[i].Y, e->gInfo[i].Z);
		ang = VerticalAng(&xyz1, xyz);
		e->gInfo[i].angH = ang;
		e->gObs[i].angH = ang;
		if (ang < (S / rou))
		{
			nGtemp -= 1;
			e->gInfo[i].flag = false;
			GPSInfo gITemp = e->gInfo[i];
			ObsData gOtemp = e->gObs[i];
			for (int r = i; r < nGtemp + 1; r++)
			{
				e->gObs[r] = e->gObs[r + 1];
				e->gInfo[r] = e->gInfo[r + 1];
			}
			e->gObs[nGtemp + 1] = gOtemp;
			e->gInfo[nGtemp + 1] = gITemp;

			temp++; i--;
		}
		else continue;
	}
	temp = 0; e->nGvalid = nGtemp;
	for (int i = 1; i < nBtemp + 1; i++)
	{
		unsigned short PRN = e->bInfo[i].PRN;
		XYZ xyz1(e->bInfo[i].X, e->bInfo[i].Y, e->bInfo[i].Z);
		ang = VerticalAng(&xyz1, xyz);
		e->bInfo[i].angH = ang;
		e->bObs[i].angH = ang;
		if (ang < (S / rou))
		{
			nBtemp -= 1;
			e->bInfo[i].flag = false;
			BDSInfo bITemp = e->bInfo[i];
			ObsData bOtemp = e->bObs[i];
			for (int r = i; r < nBtemp + 1; r++)
			{
				e->bObs[r] = e->bObs[r + 1];
				e->bInfo[r] = e->bInfo[r + 1];
			}
			e->bObs[nBtemp + 1] = bOtemp;
			e->bInfo[nBtemp + 1] = bITemp;
			temp++; i--;
		}
		else continue;
	}
	e->nBvalid = nBtemp;
	return 1;
}

/****************************
验后粗差估计排除，这里排除V>3σ的部分
对GPS信息结构体重新组织再次计算
*****************************/
int PsoteriorOutlier(ONEpoch* e, CMatrix V, double Sig0)
{
	int nGtemp = e->nGvalid; int nBtemp = e->nBvalid;
	for (int i = 0; i < V.length; i++)
	{
		if (V.value[i] > 3 * Sig0 && e->nGvalid + e->nBvalid > 5)
		{
			if (i < e->nGvalid)//GPS
			{
				nGtemp -= 1;
				e->gInfo[i].flag = false;
				GPSInfo gITemp = e->gInfo[i];
				ObsData gOtemp = e->gObs[i];
				for (int r = i; r < nGtemp + 1; r++)
				{
					e->gObs[r] = e->gObs[r + 1];
					e->gInfo[r] = e->gInfo[r + 1];
				}
				i--;
				e->gObs[nGtemp + 1] = gOtemp;
				e->gInfo[nGtemp + 1] = gITemp;
			}
			else//BDS
			{
				int itemp = i - e->nGvalid;
				nBtemp -= 1;
				e->bInfo[i].flag = false;
				BDSInfo bITemp = e->bInfo[i];
				ObsData bOtemp = e->bObs[i];
				for (int r = itemp; r < nBtemp + 1; r++)
				{
					e->bObs[r] = e->bObs[r + 1];
					e->bInfo[r] = e->bInfo[r + 1];
				}
				i--;
				e->bObs[nBtemp + 1] = bOtemp;
				e->bInfo[nBtemp + 1] = bITemp;
				//temp++;
			}
		}
		else continue;
	}
	e->nGvalid = nGtemp; e->nBvalid = nBtemp;
	return 1;
}