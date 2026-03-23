#include"decoding.h"
#include"RTK.h"
#include"lambda.h"
#include<cmath>
#include<fstream>
#include<iostream>
#include<vector>
using namespace std;
namespace {
	CMatrix Expand2FreqGeometry(const CMatrix& A)
	{
		CMatrix B = CombineRow(A, A);
		return CombineRow(B, B);
	}

	CMatrix BuildAmbiguityBlock(int nDiff, double wl1, double wl2)
	{
		CMatrix W1 = Eye(nDiff) * wl1;
		CMatrix W2 = Eye(nDiff) * wl2;
		CMatrix O = Zeros(nDiff, nDiff);
		CMatrix B1L = CombineRow(W1, O);
		CMatrix B2L = CombineRow(W2, O);
		return diag(B1L, B2L);
	}

	CMatrix BuildWeightMatrix(RTKObs* C1, RTKObs* C2, int num1, int num2)
	{
		CMatrix QL1 = GetQMatrix(C1, num1), QP1 = QL1 * 90000;
		CMatrix QL2 = GetQMatrix(C2, num2), QP2 = QL2 * 90000;
		CMatrix Q1 = diag(diag(QL1, QP1), diag(QL1, QP1));
		CMatrix Q2 = diag(diag(QL2, QP2), diag(QL2, QP2));
		return diag(Q1, Q2).inv();
	}

	CMatrix BuildResidualSystem(const RTKEpoch* E, int num, CMatrix* LL1, CMatrix* LP1, CMatrix* LL2, CMatrix* LP2, CMatrix* Ll)
	{
		*LL1 = CMatrix(E->AdrDiff1, (num - 1));
		*LP1 = CMatrix(E->PsrDiff1, (num - 1));
		*LL2 = CMatrix(E->AdrDiff2, (num - 1));
		*LP2 = CMatrix(E->PsrDiff2, (num - 1));
		CMatrix LObs = CombineRow(CombineRow(*LL1, *LP1), CombineRow(*LL2, *LP2));
		*Ll = CMatrix(E->W, 2 * (num - 1));
		*Ll = CombineRow(*Ll, *Ll);
		return LObs + (-(*Ll));
	}
}

//获取A矩阵（Xi-X0）/ρ的函数
int GetA(XYZ Sat, XYZ SatRef, XYZ Obs, double axyz[3])
{
	double distSat = CalDist(Sat, Obs);
	double distRef = CalDist(SatRef, Obs);
	double ax = (Obs.x - Sat.x) / distSat - (Obs.x - SatRef.x) / distRef;
	double ay = (Obs.y - Sat.y) / distSat - (Obs.y - SatRef.y) / distRef;
	double az = (Obs.z - Sat.z) / distSat - (Obs.z - SatRef.z) / distRef;
	axyz[0] = ax; axyz[1] = ay; axyz[2] = az;
	return 1;
}
//得到B矩阵的函数
/*
参数说明：
		x1,x2分别是基准站和流动站的坐标
		C是RTK观测值数组
		E是存储双差观测值的结构体
		wavelength是用来计算频点波长的
返回值：观测值的系数矩阵
*/
CMatrix GetBMatrix(XYZ x1, XYZ x2, RTKObs* C, RTKEpoch* E, int num)
{
	if (E->sys == GPS)
	{
		GPSInfo SR = C[E->Refnum].GRover;
		XYZ SatRef(SR.X, SR.Y, SR.Z);
		for (int i = 0; i < num; i++)
		{
			if (i == E->Refnum)continue;
			XYZ Sat(C[i].GRover.X, C[i].GRover.Y, C[i].GRover.Z); //XYZ bSat(C[i].BRover.X, C[i].BRover.Y, C[i].BRover.Z);
			GetA(Sat, SatRef, x2, C[i].matA);
		}
	}
	else if (E->sys == BDS)
	{
		BDSInfo SR = C[E->Refnum].BRover;
		XYZ SatRef(SR.X, SR.Y, SR.Z);
		for (int i = 0; i < num; i++)
		{
			if (i == E->Refnum)continue;
			XYZ Sat(C[i].BRover.X, C[i].BRover.Y, C[i].BRover.Z); //XYZ bSat(C[i].BRover.X, C[i].BRover.Y, C[i].BRover.Z);
			GetA(Sat, SatRef, x2, C[i].matA);
		}
	}
	CMatrix A;
	if (E->Refnum == 0)
	{
		A = CMatrix(C[1].matA, 3); A = A.transpose();
		for (int i = 2; i < num; i++)
		{
			CMatrix ATemp(C[i].matA, 3); ATemp = ATemp.transpose();
			A = CombineRow(A, ATemp);
		}
	}
	else
	{
		A = CMatrix(C[0].matA, 3); A = A.transpose();
		for (int i = 1; i < num; i++)
		{
			if (i != E->Refnum)
			{
				CMatrix ATemp(C[i].matA, 3); ATemp = ATemp.transpose();
				A = CombineRow(A, ATemp);
			}
		}
	}

	//直接返回设计矩阵
	return A;
}
/*
	计算协方差阵的函数，其确定方差阵的原理为：
	认为基准站和流动站的接收同卫星的观测值是等方差的，即σ(Bsae)=σ(Rover)
	而各个卫星的方差比都来自于其高度角，即σ^2(Sat)=a^2+b^2/sin^2(E)
	因此单差的方差就是(单差)=2*σ^2
	双差矩阵需要从基准站通过方差-协方差传播律得出
	Q=H*(Q0)*H'
	其中H=
	[1 0 …  0 -1
	 0 1 …  0 -1
	 ………   -1
	 0 0 …  1 -1]
	这里单系统单频的权阵就是Q阵求逆
*/
CMatrix GetQMatrix(RTKObs* C, int num)
{
	CMatrix Q(num - 1, num - 1);
	vector<double> weight(num);
	vector<double> es(num - 1);
	for (int i = 0; i < num - 1; i++)es[i] = -1;
	int j = 0; double RefSig = 0.0;
	for (int i = 0; i < num; i++)
	{
		if (C[i].Ref == true)
		{
			RefSig = C[i].sig2;
			continue;
		}
		weight[j] = C[i].sig2;
		j++;
	}
	weight[num - 1] = RefSig;
	CMatrix Q0 = diag(weight.data(), num) * 2;			//站间单差的方差矩阵
	CMatrix H0 = Eye(num - 1);
	CMatrix Es(es.data(), num - 1);
	CMatrix H = CombineCol(H0, Es);
	//MatrixShow(Q0);
	Q = H * Q0 * H.transpose();				//运用方差-协方差传播律计算权矩阵
	return Q;
}



/*
???????RTK??С??????????????????????????
?????????
	xBase??xRover-?????????????????
	C1,C2-GPS,BDS????????????
	E1,E2-????λ??????????????????????
	num1,num2-??ЧBDS,GPS????????????
	WaveLength1/2/3/4-?????????????
*/
int RTK2SysLs(XYZ xBase, XYZ* xRover, RTKObs* C1, RTKObs* C2, RTKEpoch* E1, RTKEpoch* E2,
	int num1, int num2, double WaveLength1, double WaveLength2, double WaveLength3, double WaveLength4, ofstream& F, CMatrix* DeltaX, CMatrix* Dxx, double* RMS, double* RDOP, double* Sigma0, double* Ratio)
{
	CMatrix dX(2 * num1 + 2 * num2, 1);
	CMatrix B, B1, B2, P, L, Q, A1, A2, LL1, LL11, LL2, LL22, Ll1, Ll2, LP1, LP11, LP2, LP22, LN1(num1 - 1, 1), LN2(num1 - 1, 1), LN3(num2 - 1, 1), LN4(num2 - 1, 1);
	for (int i = 0; i < num1 + 2; i++)
	{
		if (i < 3)dX.value[i] = 100;
		else dX.value[i] = 0;
	}
	while (CalNorm3(dX) >= 1e-5)
	{
		CalEpoch(GPS, xBase, *xRover, C1, E1, num1);
		CalEpoch(BDS, xBase, *xRover, C2, E2, num2);

		A1 = GetBMatrix(xBase, *xRover, C1, E1, num1);
		B1 = Expand2FreqGeometry(A1);
		CMatrix Bb1L = BuildAmbiguityBlock(num1 - 1, WaveLength1, WaveLength2);
		//B1 = CombineCol(B1,BBL);
		//???????
		//[B|??(f1)|O
		// B|O	|??(f2)]

		A2 = GetBMatrix(xBase, *xRover, C2, E2, num2);
		B2 = Expand2FreqGeometry(A2);
		CMatrix Bb2L = BuildAmbiguityBlock(num2 - 1, WaveLength3, WaveLength4);

		CMatrix O1 = Zeros(Bb1L.row, Bb2L.col);
		CMatrix O2 = Zeros(Bb2L.row, Bb1L.col);
		CMatrix BL1 = CombineCol(Bb1L, O1), BL2 = CombineCol(O2, Bb2L);
		//???????????
		//[B1,??1,O
		// B2,O,??2]?????????
		B = MatPartCombine(B1, BL1, B2, BL2);
		//MatrixShow(B);
		//MatShowFile(B,F);
		//??????????
		P = BuildWeightMatrix(C1, C2, num1, num2);
		//MatrixShow(P);

		//????в????L(L=O-C)
		CMatrix L1 = BuildResidualSystem(E1, num1, &LL1, &LP1, &LL11, &LP11, &Ll1);
		CMatrix L2 = BuildResidualSystem(E2, num2, &LL2, &LP2, &LL22, &LP22, &Ll2);
		L = CombineRow(L1, L2);
		//MatrixShow(L);

		dX = (B.transpose() * P * B).inv() * B.transpose() * P * L;
		XYZ dx(dX.value[0], dX.value[1], dX.value[2]);
		*xRover = XYZAdd(dx, *xRover);
	}
	//????????
	//dX = (B.transpose()*P*B).inv()*B.transpose()*P*L;
	CMatrix V = B * dX + (-L);
	CMatrix Qxx = (B.transpose() * P * B).inv();
	double Sig0f = ((V.transpose() * P * V) * (1.0 / (L.row - dX.row))).value[0];//????λ?????
	*RMS = sqrt((V.transpose() * V).value[0] / V.row);
	*RDOP = sqrt(tr(Get33(Qxx)));
	CMatrix Dtemp = Qxx * Sig0f;
	*Dxx = Get33(Dtemp);
	*Sigma0 = Sig0f;
	CMatrix Dx = Get33(Qxx) * Sig0f;//????????????
	//MatrixShow(Dx);
	double Sig0F = 9999;

	//lambda???й????????--??ò??????
	CMatrix Llam(dX.row - 3, 1);
	for (int i = 3; i < dX.row; i++) Llam.value[i - 3] = dX.value[i];


	CMatrix Qlam = Zeros(Qxx.row - 3, Qxx.col - 3);
	int r = 0;
	for (int i = 4; i <= Qxx.row; i++)
	{
		for (int j = 4; j <= Qxx.col; j++)
		{
			Qlam.value[r] = Qxx.find(i, j);
			r++;
		}
	}

	CMatrix Flam(2, Llam.row), slam(1, 2);
	//???????lambda????
	int info = lambda(Llam.row, 2, Llam.value, Qlam.value, Flam.value, slam.value);
	//cout << info << endl;
	//CMatrix Flam(Fl, 2 * num1 + num2 - 3, 2), slam(Sl,1,2);
	//for (int i = 0; i < dX.length-3; i++)dX.value[i + 3] = Flam.find(i+1,1);
	//MatrixShow(Flam);
	*Ratio = slam.value[1] / slam.value[0];
	if (*Ratio < 2)info = -1;
	if (info == 0)
	{
		CMatrix N(Flam.col, 1);
		for (int i = 0; i < N.row; i++)
		{
			N.value[i] = Flam.find(1, i + 1);
			if (abs(N.value[i]) < 1e-10)N.value[i] = 0;
		}
		//MatrixShow(N);
		for (int i = 0; i < num1 - 1; i++)
		{
			LN1.value[i] = LL1.value[i] - WaveLength1 * N.value[i];
			LN2.value[i] = LL11.value[i] - WaveLength2 * N.value[i + num1 - 1];
		}
		for (int i = 0; i < num2 - 1; i++)
		{
			LN3.value[i] = LL2.value[i] - WaveLength3 * N.value[i + 2 * num1 - 2];
			LN4.value[i] = LL22.value[i] - WaveLength4 * N.value[i + 2 * num1 - 2 + num2 - 1];
		}
		CMatrix LNX0 = CombineRow(LN1, LP1);
		CMatrix LNX1 = CombineRow(LN2, LP11);
		CMatrix LNX2 = CombineRow(LN3, LP2);
		CMatrix LNX3 = CombineRow(LN4, LP22);
		CMatrix LNG1 = CombineRow(LNX0, LNX1);
		CMatrix LNG2 = CombineRow(LNX2, LNX3);
 		CMatrix L0 = CombineRow(LNG1, LNG2);

		//CMatrix BN = CombineRow(B1,B2);
		double dxn[3] = { 100,100,100 };
		CMatrix dXN(dxn, 3, 1);

		while (CalNorm3(dXN) > 1e-5)
		{
			CalEpoch(GPS, xBase, *xRover, C1, E1, num1);
			CalEpoch(BDS, xBase, *xRover, C2, E2, num2);
			CMatrix BB1 = GetBMatrix(xBase, *xRover, C1, E1, num1);
			BB1 = Expand2FreqGeometry(BB1);
			CMatrix BB2 = GetBMatrix(xBase, *xRover, C2, E2, num2);
			BB2 = Expand2FreqGeometry(BB2);
			CMatrix BN = CombineRow(BB1, BB2);
			Ll1 = CMatrix(E1->W, 2 * (num1 - 1));
			Ll2 = CMatrix(E2->W, 2 * (num2 - 1));
			Ll1 = CombineRow(Ll1, Ll1);
			Ll2 = CombineRow(Ll2, Ll2);
			CMatrix Ll = CombineRow(Ll1, Ll2);
			L = L0 + (-Ll);
			//MatrixShow(BN);
			dXN = (BN.transpose() * P * BN).inv() * BN.transpose() * P * L;
			XYZ add(dXN.value[0], dXN.value[1], dXN.value[2]);
			*xRover = XYZAdd(*xRover, add);
			V = BN * dXN + (-L);
			Qxx = (BN.transpose() * P * BN).inv();
			//MatrixShow(Qxx);
			Sig0F = ((V.transpose() * P * V) * (1.0 / (L.row - dXN.row))).value[0];
		}
		*Dxx = Qxx * Sig0F;
		MatShowFile(Flam, F);
		*RMS = sqrt((V.transpose() * V).value[0] / V.row);
		*Sigma0 = Sig0F;

	}
	if (dX.value[0] > 1000)
		cout << "LS ERROR!" << endl;
	*DeltaX = dX;
	MatShowFile(dX, F);
	return info;
}

