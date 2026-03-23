#pragma once
#include<stdio.h>
#include<windows.h>

#pragma comment(lib,"WS2_32.lib")
#pragma warning(disable:4996)

struct CfgInfo
{
public:
	SOCKET Net;
	const char* NetIP = "47.114.134.129";
	int NetPort = 7190;
};

bool OpenSocket(SOCKET& sock, const char IP[], const unsigned short Port);

void CloseSocket(SOCKET& sock);

