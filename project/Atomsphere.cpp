#include"Atomsphere.h"
double HopeField(Atomsphere A, NEU neu)
{
	double E = atan(neu.u / sqrt(neu.e * neu.e + neu.n * neu.n)) * rou;
	if (E < 0)E = 0;
	double RH = St.RH * exp(-0.0006396 * (A.H - St.H));
	double  p = St.p * pow((1 - 0.0000226 * (A.H - St.H)), 5.225);
	double T = St.T - 0.0065 * (A.H - St.H);
	double e = RH * exp(-37.2465 + 0.213166 * A.T - 0.000256908 * A.T * A.T);
	double hw = 11000;
	double hd = 40136 + 148.72 * (St.T - 273.16);
	double Kw = 155.2 * 1e-7 * 4180 / (A.T * A.T) * e * (hw - A.H);
	double Kd = 155.2 * 1e-7 * p / T * (hd - A.H);
	double dTrop = Kd / sin(sqrt(E * E + 6.25) / rou) + Kw / sin(sqrt(E * E + 2.25) / rou);
	return dTrop;
}