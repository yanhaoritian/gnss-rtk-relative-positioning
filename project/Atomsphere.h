#pragma once
#include<cmath>
#include"decoding.h"

struct Atomsphere
{
	double H = 0;
	double T = 0;//温度要是绝对温标
	double p = 0;
	double RH = 0;
}const St = { 0,15 + 273.16,1013.25,0.50 };

double HopeField(Atomsphere A, NEU neu);

