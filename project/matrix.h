#pragma once
#include<cmath>
#include<iostream>
#include<stdio.h>
#include<malloc.h>
using namespace std;


/*
	矩阵类CMatrix：
	成员：
	row-行数
	col-列数

	允许的构造函数类型：
	无参数-不进行操作
	参数是一维数组-创建列向量
	参数是一维数组以及两个整形数据m,n-创建m行n列的矩阵

	成员函数
	inv()-求逆
	transpose()-求转置
	find(m,n)-找到(m,n)位置的元素
	MatrixShow(A),按照行和列显示A矩阵


	运算符的重载:
	+矩阵加法
	-矩阵乘以-1
	*乘号后跟矩阵就是矩阵乘法，跟数字就是数乘
*/

class CMatrix
{
public:
	int row = 0;
	int col = 0;
	int length = 0;
	double* value = nullptr;
	CMatrix();
	CMatrix(const CMatrix& A);
	CMatrix(int r, int c);
	CMatrix(const double* a, int l);
	CMatrix(const double* a, int row, int col);
	~CMatrix();
	double find(int rr, int cc) const;
	CMatrix transpose() const;
	CMatrix inv() const;
	CMatrix operator+(const CMatrix& B) const;
	CMatrix operator-() const;
	CMatrix operator*(const CMatrix& B) const;
	CMatrix operator*(const double& n) const;
	CMatrix& operator=(const CMatrix& B);
};

CMatrix NumMulti(const CMatrix& A, double n);//矩阵数乘函数
CMatrix Add(const CMatrix& A, const CMatrix& B);//矩阵加法函数
CMatrix MatMulti(const CMatrix& A, const CMatrix& B);//矩阵乘法函数
CMatrix Eye(int n);//创建单位阵
CMatrix CombineRow(const CMatrix& A, const CMatrix& B);//矩阵按行合并（即列数相同，B矩阵加到A后面）
CMatrix CombineCol(const CMatrix& A, const CMatrix& B);//矩阵按列合并(行数相同，B加到A后面，即C=[A,B])
void MatrixShow(const CMatrix& T);//矩阵的显示函数
CMatrix MatrixTrans(const CMatrix& A);//矩阵转置
CMatrix MatDeleteRow(const CMatrix& B, int r);//删除掉矩阵的某一行
CMatrix MatPartCombine(const CMatrix& A, const CMatrix& B, const CMatrix& C, const CMatrix& D);//矩阵分块合并
CMatrix Zeros(int r, int l);//零矩阵创建
CMatrix diag(const double* value, int n);//对角矩阵创建
CMatrix diag(const CMatrix& A, const CMatrix& B);//分块对角矩阵创建函数
void MatShowFile(const CMatrix& C, ofstream& F);//将矩阵输出到文件
CMatrix Get33(const CMatrix& A);//提取前三行前三列
double tr(const CMatrix& A);//求矩阵的迹