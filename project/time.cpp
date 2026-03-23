#include"time.h"
namespace {
	constexpr double kWeekSec = 604800.0;
	constexpr int kGpsBdsWeekOffset = 1356;
	constexpr double kGpsBdsSecOffset = 14.0;

	double NavDiffSec(const NavTime& a, const NavTime& b)
	{
		const double dweek = static_cast<double>(a.Week) - static_cast<double>(b.Week);
		return dweek * kWeekSec + (a.WeekSec - b.WeekSec);
	}

	NavTime NormalizeNavTime(int week, double sec)
	{
		while (sec >= kWeekSec) { sec -= kWeekSec; ++week; }
		while (sec < 0.0) { sec += kWeekSec; --week; }
		if (week < 0) return NavTime(0, 0.0);
		return NavTime(static_cast<unsigned short>(week), sec);
	}
}
//MJD三个构造函数，可以以三种形式来定义它，如果无参数默认为0.0
MJDTime::MJDTime()
{
	Days = 0;
	FracDay = 0.0;
	DayValue = 0.0;
}

MJDTime::MJDTime(double day)
{
	DayValue = day;
	Days = int(day);
	FracDay = day - int(day);
}

MJDTime::MJDTime(int n, double m)
{
	DayValue = n + m;
	Days = n;
	FracDay = m;
}

MJDTime UTC2MJD(const UTCTime& Utc)
{
	unsigned short year, month;
	if (Utc.Month > 2)
	{
		year = Utc.Year;
		month = Utc.Month;
	}
	else
	{
		year = Utc.Year - 1;
		month = Utc.Month + 12;
	}
	double UT = (Utc.Hour + Utc.Minute / 60.0 + Utc.Second / 3600.0) / 24.0;
	double Mjd = int(365.25 * year) + int(30.6001 * (month + 1)) + Utc.Day + UT + 1720981.5 - 2400000.5;
	MJDTime MJD(Mjd);
	return MJD;
}


/*
	简化儒略日转换为协调世界时的函数
	参数：MJD结构体
	返回值：UTC对象
*/
UTCTime MJD2UTC(const MJDTime& Mjd)
{
	double Jd = Mjd.DayValue + 2400000.5;//简化儒略日化为儒略日
	int b = int(Jd + 0.5) + 1537;
	int c = int((b - 122.1) / 365.25);
	int d = int(365.25 * c);
	int e = int((b - d) / 30.6001);
	unsigned short D = b - d - int(30.6001 * e);
	unsigned short M = e - 1 - 12 * int(e / 14.0);
	unsigned short Y = c - 4715 - int((7 + M) / 10.0);
	double day = Mjd.FracDay;
	unsigned short H = int(day * 24);
	unsigned short Min = int((day * 24 - H) * 60);
	double Sec = day * 86400 - H * 3600 - Min * 60;
	UTCTime Utc(Y, M, D, H, Min, Sec);
	return Utc;
}

/*
	简化儒略日转换为GPS时
	参数：简化儒略日结构体
	返回值：Nav时间结构体
*/
NavTime MJD2GPST(const MJDTime& Mjd)
{
	unsigned short week = int((Mjd.DayValue - 44244) / 7.0);
	double sec = (Mjd.DayValue - 44244 - 7 * week) * 86400;
	NavTime NavT(week, sec);
	return NavT;
}

/*
	Nav时转换为简化儒略日
	参数：Nav时结构体
	返回值：简化儒略日对象
*/
MJDTime GPST2MJD(const NavTime& Gpst)
{
	double Mjdays = 44244 + Gpst.Week * 7 + Gpst.WeekSec / 86400.0;
	MJDTime Mjd(Mjdays);
	return Mjd;
}




/*
UTC转换为Nav时，只需要结合前两个函数就行
*/
NavTime UTC2NavT(const UTCTime& Utc)
{
	MJDTime Mjd = UTC2MJD(Utc);
	NavTime Gpst = MJD2GPST(Mjd);
	return Gpst;
}

/********************************
Nav时的转换模块
包含：
Navminus函数
	（同样的代码重载于-当中），两个时刻的相减
		返回值是这两个时刻的差值
符号+的重载
		用于应对加上某个特定秒数时，这个秒数可正可负
*********************************/
double Navminus(const NavTime& a, const NavTime& b)
{
	return NavDiffSec(a, b);
}

double NavTime:: operator-(const NavTime& b) const
{
	return NavDiffSec(*this, b);
}


//+的重载，用于加上特定的秒数的时间转换
NavTime NavTime:: operator+(const double& t) const
{
	return NormalizeNavTime(static_cast<int>(this->Week), this->WeekSec + t);
}


/**************************************
GPS时和BDS时的相互转换，这点在双系统数据里面是很重要的
输入和输出都是NavTime类的对象
***************************************/
NavTime GPS2BDS(const NavTime& Gpst)
{
	return NormalizeNavTime(static_cast<int>(Gpst.Week) - kGpsBdsWeekOffset, Gpst.WeekSec - kGpsBdsSecOffset);
}

NavTime BDS2GPS(const NavTime& Bdst)
{
	return NormalizeNavTime(static_cast<int>(Bdst.Week) + kGpsBdsWeekOffset, Bdst.WeekSec + kGpsBdsSecOffset);
}