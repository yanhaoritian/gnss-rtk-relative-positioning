#include"matrix.h"
#include<iostream>
#include<iomanip>
#include<fstream>
//默认构造函数和析构函数都不做处理
CMatrix::CMatrix()
{
	row = 0;
	col = 0;
	length = 0;
	value = nullptr;
}
CMatrix::CMatrix(int r, int c)
{
	row = r;
	col = c;
	length = r * c;
	value = (length > 0) ? new double[length] : nullptr;
}
CMatrix::CMatrix(const CMatrix& A)
{
	row = A.row;
	col = A.col;
	length = row * col;
	value = new double[row * col];

	for (int i = 0; i < length; i++) value[i] = A.value[i];

}

CMatrix:: ~CMatrix()
{
	delete[] value;
	value = nullptr;
	row = 0;
	col = 0;
	length = 0;
}

/*
	矩阵的构造函数，如果没有行列进行控制，我们把它默认为列向量
	a-矩阵传递值的参数数组
*/
CMatrix::CMatrix(const double* a, int l)
{
	//如果没有输入行列信息，自动视为列向量
	*this = CMatrix(a, l, 1);
}

/*
	矩阵的标准构造函数
	a-传递矩阵的数值用的数组
	r-矩阵的行
	c-矩阵的列
	返回值：矩阵对象CMatrix
*/
CMatrix::CMatrix(const double* a, int r, int c)
{
	//根据行列数来确定矩阵
	row = r;
	col = c;
	length = row * col;
	if (length > 0)	value = new double[length];
	else {
		value = nullptr;
		return;
	}

	for (int i = 0; i < length; i++) value[i] = a[i];
}


/*
	Add函数，做矩阵加法
	参数：两个矩阵对象
	返回值：矩阵求和的结果
	报错类型：矩阵的维数不相同/维数是错误的
*/
CMatrix Add(const CMatrix& A, const CMatrix& B)
{
	CMatrix My(A.row, A.col);
	//矩阵的加法
	if (A.col == B.col && A.row == B.row)
	{
		for (int i = 0; i < A.length; i++)
		{
			My.value[i] = A.value[i] + B.value[i];
		}
		return My;
	}
	else //矩阵维数不一致的情况
	{
		cerr << "Not same length!";
		return CMatrix();
	}
}

/*
	矩阵的转置，是矩阵的成员函数，使用时直接在矩阵对象名后加入“.transpose()”即可
*/
CMatrix CMatrix::transpose() const
{
	CMatrix My(col, row);
	//矩阵转置

	int k, l, x, y;
	for (x = 0; x < row; x++)
	{
		for (y = 0; y < col; y++)
		{
			k = x * col + y;
			l = y * row + x;
			My.value[l] = value[k];
		}
	}
	//CMatrix Res(bvalue, rr, ll);
	//bvalue = nullptr;
	return My;
}


/*
	矩阵的乘法，使用*进行重载也可以达到相同的效果，参数A和B分别是待相乘的两个矩阵
*/
CMatrix MatMulti(const CMatrix& A, const CMatrix& B)
//矩阵乘法的函数，左乘右乘有顺序
{
	CMatrix Res(A.row, B.col);
	if (A.col != B.row)
		//如果出现矩阵的维数有问题，即报错
	{
		cerr << "MultipyERR: Incorrect rows/columns!" << endl;
	}
	//int resRow=A.row;	
	//int resCol = B.col;								//矩阵的行和列由相乘的两个矩阵行列决定
	//double* resValue=new double[resRow*resCol];					//结果矩阵的元素在这里储存
	int  i, j, k;
	for (i = 0; i < A.row; i++)
	{
		for (j = 0; j < B.col; j++)
		{
			Res.value[i * B.col + j] = 0;
			for (k = 0; k < A.col; k++)
			{
				Res.value[i * B.col + j] = Res.value[i * B.col + j] + A.value[i * A.col + k] * B.value[k * B.col + j];
			}
		}
	}
	//CMatrix resMat(resValue,resRow,resCol);
	//resValue = nullptr;
	return Res;
}


/*
	寻找矩阵中的某个元素
	rr-该元素所在的行
	cc-该元素所在的列
	返回值：矩阵元素A(rr,cc)
*/
double CMatrix::find(int rr, int cc) const
{
	rr = rr - 1;
	cc = cc - 1;
	if (rr > row || rr < 0 || cc > col || cc < 0)
	{
		cerr << "Incorrect element searching!" << endl;
		return 0.0;
	}
	else
	{
		int num = rr * col + cc;
		double elementFind = value[num];
		return elementFind;
	}
}

/*
	矩阵求逆函数，直接使用矩阵对象.inv()即可
	参数：无，要求必须是方阵
*/
CMatrix CMatrix::inv() const
{
	CMatrix Res(row, col);
	if (this->col != this->row || this->col == 0 || this->row == 0)
	{
		cerr << "InvERR: Not a square matrix!" << endl;
		CMatrix Err; return Err;
		//exit(1);
	}
	int n = this->row;
	double* a = this->value;
	double* b = Res.value;
	int i, j, k, l, u, v;				//均为循环变量
	int is[500], js[1500];   /* matrix dimension  (有限制)*/
	double d, p;

	/* 将输入矩阵赋值给输出矩阵b，下面对b矩阵求逆，a矩阵不变 */
	for (k = 0; k < n * n; k++)
	{
		b[k] = a[k];
	}
	for (k = 0; k < n; k++)
	{
		d = 0.0;
		for (i = k; i < n; i++)   /* 查找右下角方阵中主元素的位置 */
		{
			for (j = k; j < n; j++)
			{
				l = n * i + j;
				p = fabs(b[l]);
				if (p > d)
				{
					d = p;
					is[k] = i;
					js[k] = j;
				}
			}
		}

		if (d < 1.0E-12)
		{
			printf("Divided by 0 in MatrixInv!\n");
			CMatrix Btemp;
			return Btemp;
			//return;
			//exit(1);
		}

		if (is[k] != k)  /* 对主元素所在的行与右下角方阵的首行进行调换 */
		{
			for (j = 0; j < n; j++)
			{
				u = k * n + j;
				v = is[k] * n + j;
				p = b[u];
				b[u] = b[v];
				b[v] = p;
			}
		}

		if (js[k] != k)  /* 对主元素所在的列与右下角方阵的首列进行调换 */
		{
			for (i = 0; i < n; i++)
			{
				u = i * n + k;
				v = i * n + js[k];
				p = b[u];
				b[u] = b[v];
				b[v] = p;
			}
		}

		l = k * n + k;
		b[l] = 1.0 / b[l];  /* 初等行变换 */
		for (j = 0; j < n; j++)
		{
			if (j != k)
			{
				u = k * n + j;
				b[u] = b[u] * b[l];
			}
		}
		for (i = 0; i < n; i++)
		{
			if (i != k)
			{
				for (j = 0; j < n; j++)
				{
					if (j != k)
					{
						u = i * n + j;
						b[u] = b[u] - b[i * n + k] * b[k * n + j];
					}
				}
			}
		}
		for (i = 0; i < n; i++)
		{
			if (i != k)
			{
				u = i * n + k;
				b[u] = -b[u] * b[l];
			}
		}
	}

	for (k = n - 1; k >= 0; k--)  /* 将上面的行列调换重新恢复 */
	{
		if (js[k] != k)
		{
			for (j = 0; j < n; j++)
			{
				u = k * n + j;
				v = js[k] * n + j;
				p = b[u];
				b[u] = b[v];
				b[v] = p;
			}
		}
		if (is[k] != k)
		{
			for (i = 0; i < n; i++)
			{
				u = i * n + k;
				v = is[k] + i * n;
				p = b[u];
				b[u] = b[v];
				b[v] = p;
			}
		}
	}
	//CMatrix inverse(b,n,n);
	//b = nullptr;
	return Res;
}

CMatrix Eye(int n)
//构建一个n*n的单位矩阵
{
	if (n <= 0)
	{
		cerr << "EyeErr: n<=0！" << endl;
		//return;
		//exit(1);
	}
	CMatrix Eye(n, n);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (j == i)
			{
				Eye.value[i * n + j] = 1;
			}
			else
			{
				Eye.value[i * n + j] = 0;
			}
		}
	}
	//CMatrix Eye(Res,n,n);
	//Res = nullptr;
	return Eye;
}

CMatrix NumMulti(const CMatrix& A, double n)//矩阵数乘函数
{
	CMatrix Res(A.row, A.col);
	//double* temp = new double[A.length];
	for (int i = 0; i < A.length; i++)
	{
		Res.value[i] = A.value[i] * n;
	}
	return Res;
}

/*
矩阵按行合并（即要求列数相同，B矩阵加到A后面）
同matlab中的C=[A ; B]
*/
CMatrix CombineRow(const CMatrix& A, const CMatrix& B)
{
	if (A.col != B.col)
	{
		cerr << "CombineERR: Not Same Columns!" << endl;
		CMatrix nErr;
		return nErr;
		//exit(1);
	}
	CMatrix C(A.row + B.row, A.col);
	//double *cvalue = new double[A.length + B.length];


	for (int i = 0; i < A.length; i++)
	{
		C.value[i] = A.value[i];
	}
	for (int j = 0; j < B.length; j++)
	{
		C.value[A.length + j] = B.value[j];
	}
	//CMatrix C(cvalue,A.row+B.row,A.col);
	//cvalue = nullptr;
	return C;
}


//矩阵按列合并(行数相同，B加到A后面，即C=[A,B])
CMatrix CombineCol(const CMatrix& A, const CMatrix& B)
{
	if (A.row != B.row)
	{
		cerr << "CombineERR: Not Same Rows!" << endl;
		//exit(1);
	}
	CMatrix Ctemp = CombineRow(A.transpose(), B.transpose());
	CMatrix C = Ctemp.transpose();
	return C;
}

//矩阵的显示函数，将按照行列进行表示
void MatrixShow(const CMatrix& T)
{
	int num = 0;
	for (int tem = 0; tem < T.row; tem++)
	{
		for (int rem = 0; rem < T.col; rem++)
		{
			cout << T.value[num] << " ";
			num++;
		}
		cout << endl;
	}
}


/*
  函数名称:MatrixTrans
  函数功能:矩阵转置
  参数:
  A-输入矩阵A
*/
CMatrix MatrixTrans(const CMatrix& A)
{
	CMatrix Res(A.col, A.row);
	//int rr=0,ll = 0;
	//double *bvalue = new double[A.length];
	int k, l, x, y;
	for (x = 0; x < A.row; x++)
	{
		for (y = 0; y < A.col; y++)
		{
			k = x * A.col + y;
			l = y * A.row + x;
			Res.value[l] = A.value[k];
		}
	}
	//rr = A.col;
	//ll = A.row;
	//CMatrix Res(bvalue,rr,ll);
	//bvalue = nullptr;
	return Res;
}

CMatrix CMatrix::operator+(const CMatrix& B) const
{
	//矩阵加法的重载符
	if (this->col == B.col && this->row == B.row)
	{
		CMatrix Res(row, col);
		//int col = this->col;
		//int row = this->row;
		//double* c = new double[this->length];
		for (int i = 0; i < this->length; i++)
		{
			Res.value[i] = this->value[i] + B.value[i];
		}
		//CMatrix Res(c, row, col);
		//c = nullptr;
		return Res;
	}
	else //矩阵维数不一致的情况直接报错
	{
		cerr << "PlusERR: Not same length!";
		return CMatrix();
	}
}

CMatrix CMatrix::operator-() const
{
	//矩阵取负号的重载
	CMatrix Minus(row, col);
	//double* minus = new double[this->length];
	for (int i = 0; i < length; i++)
	{
		Minus.value[i] = -value[i];
	}
	//CMatrix Minus(minus,row,col);
	//minus = nullptr;
	return Minus;
}

CMatrix CMatrix::operator*(const CMatrix& B) const
{
	if (this->col != B.row)
		//如果出现矩阵的维数有问题，即报错
	{
		cerr << "MultipyERR: Incorrect rows/columns!" << endl;
		return CMatrix();
	}
	else
	{
		CMatrix Res(this->row, B.col);								//矩阵的行和列由相乘的两个矩阵行列决定
		//double* resValue = new double[resRow*resCol];					//结果矩阵的元素在这里储存
		int  i, j, k;
		for (i = 0; i < this->row; i++)
		{
			for (j = 0; j < B.col; j++)
			{
				Res.value[i * B.col + j] = 0;
				for (k = 0; k < this->col; k++)
				{
					Res.value[i * B.col + j] = Res.value[i * B.col + j] + this->value[i * this->col + k] * B.value[k * B.col + j];
				}
			}
		}
		//CMatrix resMat(resValue, resRow, resCol);
		//resValue = nullptr;
		return Res;
	}
}

CMatrix CMatrix::operator*(const double& n) const//矩阵数乘函数
{
	CMatrix Temp(row, col);
	//int m = length;
	//double* temp = new double[m];
	for (int i = 0; i < length; i++)
	{
		double a = value[i];
		Temp.value[i] = a * n;
	}
	//CMatrix Res(temp, row, col);
	//temp = nullptr;
	return Temp;
}

CMatrix& CMatrix::operator=(const CMatrix& A)
{
	if (this == &A) return *this;
	delete[] value;
	row = A.row;
	col = A.col;
	length = row * col;
	value = (length > 0) ? new double[length] : nullptr;

	for (int i = 0; i < length; i++) value[i] = A.value[i];
	return *this;
}

/*
矩阵删除某一行的函数，剔掉个别不健康观测值时使用
参数：B需要删除行的矩阵，r删除的行的标号
*/
CMatrix MatDeleteRow(const CMatrix& B, int r)
{
	if (B.row > r && B.row > 0 && B.col > 0)
	{
		CMatrix C(B.row - 1, B.col);
		for (int s = 0; s < (r - 1) * B.col; s++) C.value[s] = B.value[s];
		for (int s = (r - 1) * B.col; s < C.length; s++) C.value[s] = B.value[s + B.col];
		return C;
	}
	else
	{
		CMatrix C = B;
		return C;
	}
}

/*
矩阵合并函数，ABCD是矩阵的四个块
合并方式为：
[A|B
 C|D]
*/
CMatrix MatPartCombine(const CMatrix& A, const CMatrix& B, const CMatrix& C, const CMatrix& D)
{
	CMatrix AC = CombineRow(A, C);
	CMatrix BD = CombineRow(B, D);
	CMatrix ABCD = CombineCol(AC, BD);
	return ABCD;
}

//创建一个全零矩阵
CMatrix Zeros(int r, int l)
{
	CMatrix Zero(r, l);
	for (int i = 0; i < r * l; i++) Zero.value[i] = 0;
	return Zero;
}

//创建一个对角矩阵，其中value只存储对角值
CMatrix diag(const double* value, int n)
{
	CMatrix D(n, n);
	int r = 0;
	for (int i = 0; i < n * n; i++)
	{
		if (i == (n + 1) * r)
		{
			D.value[i] = value[r];
			r++;
		}
		else D.value[i] = 0;
	}
	return D;
}

//创建一个分块对角矩阵
//[A|O
// O|B]
CMatrix diag(const CMatrix& A, const CMatrix& B)
{
	CMatrix O1 = Zeros(A.row, B.col);
	CMatrix O2 = Zeros(B.row, A.col);
	return MatPartCombine(A, O1, O2, B);
}

//输出矩阵值到文件中的函数
void MatShowFile(const CMatrix& T, ofstream& F)
{
	int num = 0;
	for (int tem = 0; tem < T.row; tem++)
	{
		for (int rem = 0; rem < T.col; rem++)
		{
			F << setiosflags(ios::fixed) << setprecision(4) << T.value[num] << "\t";
			num++;
		}
		F << endl;
	}
}


//矩阵的对角化，D是对角化的矩阵
/* LD factorization (Q=L*diag(D)*L') -----------------------------------------*/
static int LD(const CMatrix Q, CMatrix L, CMatrix D)   //此处分解算法同 FMFAC5（参见"模糊度经典版.pdf"） 
{
	int n = Q.col;
	int i, j, k, info = 0;
	double a; CMatrix A(n, n); //mat():给矩阵分配内存
	A = Q;
	//memcpy(A, Q, sizeof(double)*n*n); //将Q阵复制到A中
	for (i = n - 1; i >= 0; i--) {          //FMFAC5中是从i=n开始的，为何有区别？由于矩阵元素是从0到n*n-1，而FMFAC5中是1到n*n，故此处需要多减1
		if ((D.value[i] = A.value[i + i * n]) <= 0.0) { info = -1; break; }  //D(i,i)=Q(i,i),注意此处D是一个列向量而非矩阵
		a = sqrt(D.value[i]);                                //a=sqrt(Q(i,i))
		for (j = 0; j <= i; j++) L.value[i + j * n] = A.value[i + j * n] / a;      //L(i,1:i)=Q(i,1:i)/sqrt(Q(i,i))
		for (j = 0; j <= i - 1; j++) for (k = 0; k <= j; k++) A.value[j + k * n] -= L.value[i + k * n] * L.value[i + j * n];  //Q(j,1:j)=Q(j,1:j)-L(i,1:j)L(i,j) 同样的由于矩阵原因，j从0开始   
		for (j = 0; j <= i; j++) L.value[i + j * n] /= L.value[i + i * n];       //L(i,1:i)=L(i,1:i)/L(i,i)
	}
	//free(A);
	L = L.transpose();
	if (info) printf("LD factorization error\n");
	return info;
}

//获取前三行前三列
CMatrix Get33(const CMatrix& A)
{
	if (A.col < 3 || A.row < 3)
	{
		cout << "MatrixErr:Not a 3*3Mat." << endl;
		return A;
	}
	CMatrix B(3, 3);
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
			B.value[3 * i + j] = A.value[3 * i + j];
	}
	return B;
}

double tr(const CMatrix& A)
{
	double tr = 0;
	if (A.col != A.row || A.length < 0)
	{
		cout << "ERROR: Mat Tr not exist!" << endl;
		return tr;
	}
	for (int i = 1; i <= A.row; i++)tr = tr + A.find(i, i);
	return tr;
}