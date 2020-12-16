#include <iostream>
#include <vector>
#include <algorithm>
#include <omp.h>
using namespace std;
void ModelGenerationCUBE(double SideLength, int CutNum, int nn, int NLeng, int Nx, int Ny, int Nz, double** gcoord, int** nodes, int wmth);
int main()
{
	int CutNum, nn, NLeng, Nx, Ny, Nz;
	CutNum = 50;       //切割次数
	double Magnitude = 1;     //周期幅度
	double BasicR = 1;        //纤维半径
	double Cycle = 10;        //纤维周期长度
	int Curve = 1;
	int wmth = 1;      //编织方式   0为1D，1为2D，2为2.5D，3为3D
	nn = CutNum + 1;
	double SideLength = 10;
	//参数化部分//
	double CrossSectExp = 2;
	double p = 1;//0.03413
	double q = 1;//0.2
   ///设定xyz的比例（晶格大小）///
	double index_xyz = 1;// 0.1;
	///////////////////////////////
	NLeng = nn;
	Nx = NLeng;
	Ny = NLeng;
	Nz = NLeng;

	switch (wmth)
	{
		case(0):
			Nz = 0.5 * (NLeng - 1) + 1;
			break;
		case(1):
			Nz = 0.5 * (NLeng - 1) + 1;
			break;
		case(2):
			Magnitude = 0.9;     //周期幅度
			BasicR = 0.2;        //纤维半径
		    Cycle = 6.72;        //纤维周期长度
			SideLength = 10 * 0.672;
			Ny = ceil(NLeng * (3.6 / 6.72)) + 1;
			Nz = ceil(NLeng * (4.05 / 6.72)) + 1;
			break;
		case(3):
			Nz = NLeng;
			break;
	}
	//////////////////////////////////////////////
	//vector< vector<double> > gcoord((Nx * Ny * Nz), vector<double>(3, 0));
	//vector< vector<int> > nodes(((Nx - 1) * (Ny - 1) * (Nz - 1) * 5), vector<int>(4, 0));
	double** gcoord = new double* [(Nx * Ny * Nz)];
	for (int i = 0; i < (Nx * Ny * Nz); i++)
		gcoord[i] = new double[3]();
	int** nodes = new int* [((Nx - 1) * (Ny - 1) * (Nz - 1) * 5)];
	for (int i = 0; i < ((Nx - 1) * (Ny - 1) * (Nz - 1) * 5); i++)
		nodes[i] = new int[4]();
	ModelGenerationCUBE(SideLength, CutNum, nn, NLeng, Nx, Ny, Nz, gcoord, nodes, wmth);
	int IEn, IE2n, FEn;
	switch (wmth)
	{
	case(0):
		IEn = 3, IE2n = 3, FEn = 6;
		break;
	case(1):
		IEn = 3, IE2n = 3, FEn = 6;
		break;
	case(2):
		IEn = 12, IE2n = 12, FEn = 24;
		break;
	case(3):
		IEn = 8, IE2n = 1, FEn = 8;
		break;
	}

	int** InnerElementsmain = new int* [IEn];
	for (int i = 0; i < IEn; i++)
		InnerElementsmain[i] = new int[1]();

	int** InnerElementsmain2st = new int* [IE2n];
	for (int i = 0; i < IE2n; i++)
		InnerElementsmain2st[i] = new int[1]();

	int* OuterElementsmain = new int[1];

	int** FiberEle = new int* [FEn];
	for (int i = 0; i < FEn; i++)
		FiberEle[i] = new int[1]();







	for (int i = 0; i < (Nx * Ny * Nz); i++)
		delete[]gcoord[i];
	delete[]gcoord;
	//
	for (int i = 0; i < ((Nx - 1) * (Ny - 1) * (Nz - 1) * 5); i++)
		delete[]nodes[i];
	delete[]nodes;
	//
	for (int i = 0; i < IEn; i++)
		delete[]InnerElementsmain[i];
	delete[]InnerElementsmain;
	//
	for (int i = 0; i < IE2n; i++)
		delete[]InnerElementsmain2st[i];
	delete[]InnerElementsmain2st;
	//
	delete[]OuterElementsmain;
	//
	for (int i = 0; i < FEn; i++)
		delete[]FiberEle[i];
	delete[]FiberEle;
	int aa = 55;
}




