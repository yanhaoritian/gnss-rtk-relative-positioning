// project.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
#include<iostream>
#include"fstream"
#include<iomanip>
#include"Atomsphere.h"
#include"coordinate.h"
#include"Matrix.h"
#include"decoding.h"
#include"RTK.h"
#pragma warning(disable : 4996)
using namespace std;

ONEpoch EBase, ERover;

int main()
{
	
	//FILE* FRover = fopen("oem719-202203031500-1.bin", "rb");//零基线
	//FILE* FBase = fopen("oem719-202203031500-2.bin", "rb");

	//FILE* FRover = fopen("oem719-202203170900-1.bin", "rb");//短基线
	//FILE* FBase = fopen("oem719-202203170900-2.bin", "rb");

	
	unsigned char buffBase[MAXRAWLEN] = { 0 };					//基准站数组
	unsigned char buffRover[MAXRAWLEN] = { 0 };					//流动站数组
	CfgInfo Sect1, Sect2;
	Sect1.NetIP = "47.114.134.129"; Sect2.NetIP = "8.140.46.126";
	Sect1.NetPort = 7190; Sect2.NetPort = 3002;
	RTKSockets(Sect1, Sect2);
	//PPKFile(FBase, FRover, buffBase, buffRover, &EBase, &ERover);
	return 0;
}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
