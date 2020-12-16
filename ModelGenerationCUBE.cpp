#include <iostream>
//#include <vector>
//#include <omp.h>
using namespace std;
void ModelGenerationCUBE(double SideLength, int CutNum, int nn, int NLeng, int Nx, int Ny, int Nz, double** gcoord, int** nodes, int wmth)
{
	double  Lx, Ly, Lz, hx, hy, hz;

	Lx = 1; //Nx = 10 * NLeng - 9;
	Ly = 1; //Ny = 2 * NLeng - 1;
	Lz = 1; //Nz = NLeng;

	double zy;
	zy = Lx / 9;
	hx = Lx / (Nx - 1);
	//vector<double> x(Nx, 0);
	int* x = new int[Nx];
	//omp_set_num_threads(8);
	//#pragma omp parallel for
	for (int i = 0; i < (Nx - 1); i++)
	{
		x[i + 1] = x[i] + hx;
	}

	hy = Ly / (Ny - 1);
	//vector<double> y(Ny, 0);
	int* y = new int[Ny];
	//#pragma omp parallel for
	for (int i1 = 0; i1 < (Ny - 1); i1++)
	{
		y[i1 + 1] = y[i1] + hy;
	}

	hz = Lz / (Nz - 1);
	//vector<double> z(Nz, 0);
	int* z = new int[Nz];
	//#pragma omp parallel for
	for (int i2 = 0; i2 < (Nz - 1); i2++)
	{
		z[i2 + 1] = z[i2] + hz;
	}

#pragma omp parallel for
	for (int i = 0; i < (Ny * Nz); i++)
	{
		for (int j = 0; j < Nx; j++)
		{
			gcoord[j + i * Nx][0] = x[j];
		}
	}
	//
	for (int i = 0; i < Nz; i++)
	{
		for (int i2 = 0; i2 < Ny; i2++)
		{
			for (int j = 0; j < Nx; j++)
			{
				gcoord[j + i2 * Nx + i * Nx * Ny][1] = y[i2];
			}
		}
	}
#pragma omp parallel for
	for (int i = 0; i < Nz; i++)
	{
		for (int j = 0; j < (Nx * Ny); j++)
		{
			gcoord[j + i * Nx * Ny][2] = z[i];
		}
	}

	//vector<vector<int>> type1(Nx - 1, vector<int>(Ny - 1, 0));
	int** type1 = new int* [(Nx - 1)];
	for (int i = 0; i < (Nx - 1); i++)
		type1[i] = new int[Ny - 1];
	//#pragma omp parallel for
	for (int i = 0; i < (Nx - 1); i = i + 2)
	{
		for (int j = 0; j < (Ny - 1); j = j + 2)
		{
			type1[i][j] = 1;
		}
	}
	//#pragma omp parallel for
	for (int i = 1; i < (Nx - 1); i = i + 2)
	{
		for (int j = 1; j < (Ny - 1); j = j + 2)
		{
			type1[i][j] = 1;
		}
	}


	//vector<vector<int>> type2(Nx - 1, vector<int>(Ny - 1, 0));
	int** type2 = new int* [(Nx - 1)];
	for (int i = 0; i < (Nx - 1); i++)
		type2[i] = new int[Ny - 1];
	//#pragma omp parallel for
	for (int i = 0; i < (Nx - 1); i = i + 2)
	{
		for (int j = 1; j < (Ny - 1); j = j + 2)
		{
			type2[i][j] = 1;
		}
	}
	//#pragma omp parallel for
	for (int i = 1; i < (Nx - 1); i = i + 2)
	{
		for (int j = 0; j < (Ny - 1); j = j + 2)
		{
			type2[i][j] = 1;
		}
	}

	//vector <vector<vector<int>>> Type(Nz - 1, vector<vector<int>>(Nx - 1, vector<int>(Ny - 1, 0)));
	int*** Type = new int** [Nz - 1];
	for (int i = 0; i < (Nz - 1); i++)
	{
		Type[i] = new int* [Nx - 1];
	}
	for (int i = 0; i < (Nz - 1); i++)
	{
		for (int j = 0; j < (Nx - 1); j++)
		{
			Type[i][j] = new int[Ny - 1];
		}
	}



#pragma omp parallel for
	for (int i = 0; i < (Nz - 1); i++)
	{
		if ((i % 2) == 0)
		{
			Type[i] = type1;
		}
		else
		{
			Type[i] = type2;
		}
	}

	//vector<vector<int>> Tab_connect((Nx - 1) * (Ny - 1) * (Nz - 1) * 5, vector<int>(4, 0));
	int Ne = 0;
	int pts[8] = { 0 };
	//#pragma omp parallel for
	for (int k = 1; k < Nz; k++)
	{
		for (int j = 1; j < Ny; j++)
		{
			for (int i = 1; i < Nx; i++)
			{
				//pts=(k-1)*Nx*Ny+(j-1)*Nx+[i i+1 Nx+i+1 Nx+i Nx* Ny + [i i + 1 Nx + i + 1 Nx + i]]
				pts[0] = (k - 1) * Nx * Ny + (j - 1) * Nx + i - 1;
				pts[1] = (k - 1) * Nx * Ny + (j - 1) * Nx + i + 1 - 1;
				pts[2] = (k - 1) * Nx * Ny + (j - 1) * Nx + Nx + i + 1 - 1;
				pts[3] = (k - 1) * Nx * Ny + (j - 1) * Nx + Nx + i - 1;
				pts[4] = (k - 1) * Nx * Ny + (j - 1) * Nx + Nx * Ny + i - 1;
				pts[5] = (k - 1) * Nx * Ny + (j - 1) * Nx + Nx * Ny + i + 1 - 1;
				pts[6] = (k - 1) * Nx * Ny + (j - 1) * Nx + Nx * Ny + Nx + i + 1 - 1;
				pts[7] = (k - 1) * Nx * Ny + (j - 1) * Nx + Nx * Ny + Nx + i - 1;

				if (Type[k - 1][i - 1][j - 1] == 1)
				{
					nodes[Ne][0] = pts[0], nodes[Ne][1] = pts[1], nodes[Ne][2] = pts[2], nodes[Ne][3] = pts[5];
					nodes[Ne + 1][0] = pts[2], nodes[Ne + 1][1] = pts[3], nodes[Ne + 1][2] = pts[0], nodes[Ne + 1][3] = pts[7];
					nodes[Ne + 2][0] = pts[2], nodes[Ne + 2][1] = pts[7], nodes[Ne + 2][2] = pts[5], nodes[Ne + 2][3] = pts[6];
					nodes[Ne + 3][0] = pts[0], nodes[Ne + 3][1] = pts[5], nodes[Ne + 3][2] = pts[7], nodes[Ne + 3][3] = pts[4];
					nodes[Ne + 4][0] = pts[0], nodes[Ne + 4][1] = pts[2], nodes[Ne + 4][2] = pts[7], nodes[Ne + 4][3] = pts[5];
				}
				else
				{
					nodes[Ne][0] = pts[0], nodes[Ne][1] = pts[1], nodes[Ne][2] = pts[3], nodes[Ne][3] = pts[4];
					nodes[Ne + 1][0] = pts[1], nodes[Ne + 1][1] = pts[2], nodes[Ne + 1][2] = pts[3], nodes[Ne + 1][3] = pts[6];
					nodes[Ne + 2][0] = pts[1], nodes[Ne + 2][1] = pts[6], nodes[Ne + 2][2] = pts[4], nodes[Ne + 2][3] = pts[5];
					nodes[Ne + 3][0] = pts[3], nodes[Ne + 3][1] = pts[4], nodes[Ne + 3][2] = pts[6], nodes[Ne + 3][3] = pts[7];
					nodes[Ne + 4][0] = pts[3], nodes[Ne + 4][1] = pts[1], nodes[Ne + 4][2] = pts[6], nodes[Ne + 4][3] = pts[4];
				}
				Ne = Ne + 5;
			}
		}
	}

	//#pragma omp parallel for
	//	for (int m = 0; m < nodes.size(); m++)
	//	{
	//		for (int n = 0; n < nodes[m].size(); n++)
	//		{
	//			nodes[m][n] = Tab_connect[m][n]-1;
	//		}
	//	}
#pragma omp parallel for
	for (int i = 0; i < Nx * Ny * Nz; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			gcoord[i][j] = gcoord[i][j] * SideLength - 0.5 * SideLength;
		}
		cout << endl;
	}

	//#pragma omp parallel for
	//	for (int i = 0; i < Nx * Ny * Nz; i++)
	//	{
	//			gcoord[i][2] = gcoord[i][2]* 0.5 ;
	//	}
	if (wmth == 0 || wmth == 1)
	{
#pragma omp parallel for
		for (int i = 0; i < Nx * Ny * Nz; i++)
		{
			gcoord[i][2] = gcoord[i][2] * 0.5;
		}
	}
	if (wmth == 2)
	{
#pragma omp parallel for
		for (int i = 0; i < Nx * Ny * Nz; i++)
		{
			gcoord[i][1] = gcoord[i][1] * 3.6 / 6.72;
			gcoord[i][2] = gcoord[i][2] * 4.05 / 6.72;
		}
	}
	int aa = 55;
	delete[]x;
	delete[]y;
	delete[]z;
	//for (int i = 0; i < (Nx * Ny * Nz); i++)
	//	delete[]coor[i];
	//delete[]coor;//删除数组
	for (int i = 0; i < (Nx - 1); i++)
		delete[]type1[i];
	delete[]type1;//删除数组
	for (int i = 0; i < (Nx - 1); i++)
		delete[]type2[i];
	delete[]type2;//删除数组

	delete[]Type;
	//for (int x = 0; x < (Nz-1); x++)
	//{
	//	for (int y = 0; y < (Nx - 1); y++)
	//	{
	//		delete[] Type[x][y];//释放Z这一层
	//	}
	//}
	//for (int x = 0; x < (Nz - 1); x++)
	//{
	//	delete[] Type[x];//释放Y这一层
	//}
	//delete[] Type;//释放X这一层
	//for (int i = 0; i < ((Nx - 1) * (Ny - 1) * (Nz - 1) * 5); i++)
	//	delete[]Tab_connect[i];
	//delete[]Tab_connect;//删除数组

}