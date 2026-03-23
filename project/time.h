#pragma once
#include<iostream>
#include<cmath>
using namespace std;

// 简化儒略日容器。
// Days：整数日，FracDay：小数日，DayValue：总日数。
class MJDTime
{
public:
	int Days;
	double FracDay;
	double DayValue;

	MJDTime();
	MJDTime(double day);
	MJDTime(int n, double m);
};

// 导航时间容器（周 + 周内秒）。
class NavTime
{
public:
	unsigned short Week;
	double WeekSec;

	NavTime()
	{
		Week = 0;
		WeekSec = 0.0;
	}
	NavTime(unsigned short week, double sec)
	{
		Week = week;
		WeekSec = sec;
	}
	double operator-(const NavTime& b) const;
	NavTime operator+(const double& t) const;
};

// UTC 时间容器。
class UTCTime
{
public:
	unsigned short Year;
	unsigned short Month;
	unsigned short Day;
	unsigned short Hour;
	unsigned short Minute;
	double Second;

	UTCTime() { Year = 0; Month = 0; Day = 0; Hour = 0; Minute = 0; Second = 0.0; }

	UTCTime(unsigned short yy, unsigned short mm, unsigned short dd, unsigned short hh, unsigned short min, double s)
	{
		Year = yy; Month = mm; Day = dd; Hour = hh; Minute = min; Second = s;
	}
};

// 时间转换辅助函数。
MJDTime UTC2MJD(const UTCTime& Utc);
UTCTime MJD2UTC(const MJDTime& Mjd);
NavTime MJD2GPST(const MJDTime& Mjd);
MJDTime GPST2MJD(const NavTime& Gpst);
NavTime UTC2NavT(const UTCTime& Utc);
NavTime GPS2BDS(const NavTime& Gpst);
NavTime BDS2GPS(const NavTime& Bdst);
double Navminus(const NavTime& a, const NavTime& b);
