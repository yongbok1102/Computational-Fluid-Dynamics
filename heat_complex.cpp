#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
int min(int x, int y)
{
	if(x>y)
		return y;
	else
		return x;
}

/*----------------------------------------------*/
/*source code for steady-state heat transfer	*/
/*through the complex plate			*/
/*Method for solving PDE: FDM			*/
/*Method for solving linear system: Gauss-Seidel*/
/*Written by yongbok1102			*/
/*----------------------------------------------*/

int main()
{

	//grid generation
	int nx1 = 150; int ny1 = 100;
	int nx2 = 100; int ny2 = 50;

	int np1 = (1 + nx1)*(1 + ny1);
	int np2 = (1 + nx2)*ny2;
	int np = np1 + np2;
	int ncell = nx1*ny1 + nx2*ny2;
	double dx1 = 1.5 / nx1; double dy1 = 1. / ny1;
	double dx2 = 1. / nx2; double dy2 = 0.5 / ny2;

	double* X; X = new double[np];
	double* Y; Y = new double[np];
	double* T; T = new double[np];

	int n = 0;

	//numbering the points
	for (int j = 0; j <= ny1; j++)
	{
		for (int i = 0; i <= nx1; i++)
		{
			X[n] = i*dx1;
			Y[n] = j*dy1;
			n++;
		}
	}
	for (int j = 1; j <= ny2; j++)
	{
		for (int i = 0; i <= nx2; i++)
		{
			X[n] = i*dx2;
			Y[n] = 1.0 + j*dy2;
			n++;
		}
	}

	//initialization
	for (int i = 0; i < np; i++)
	{
		T[i] = 0;
	}

	//Iterative calculation
	for (int itr = 1; itr <= 100000; itr++)
	{
		cout << "iteration " << itr << endl;
		

		for (int i = 0; i <np; i++)
		{	
			//region 1
			if(i<np1)
			{
				//left side (BC)
				if(abs(X[i])<1e-006)
					T[i]=100;
				//right side (BC)
				else if(abs(X[i]-1.50)<1e-006)
					T[i]=T[i-1];
				//lower side (BC)
				else if(abs(Y[i])<1e-006)
					T[i]=100;
				//region 1(excluding the boundaries)
				else if(abs(X[i])>1e-006 && abs(X[i]-1.50)>1e-006 && abs(Y[i])>1e-006 && abs(Y[i]-1.0)>1e-006)
					T[i]=0.25*(T[i-1]+T[i+1]+T[i-nx1-1]+T[i+nx1+1]);
	
				//between region 1 and 2
				else
				{
					if(i==ny1*(1+nx1))
						T[i]=100;
					//upper side(BC)
					else if(i>=ny1*(1+nx1)+min(nx1,nx2) && i<np1)
						T[i]=200;
					//close to region 2
					else
						T[i]=0.25*(T[i-1]+T[i+1]+T[i-nx1-1]+T[i+nx1+1]);
				}
			}
			//region 2
			else
			{
				//left side(BC)
				if(abs(X[i])<1e-006)
					T[i]=100;
				//right side(BC)
				else if(abs(X[i]-1)<1e-006)
					T[i]=200;
				//upper side(BC)
				else if(abs(Y[i]-1.5)<1e-006)
					T[i]=T[i-nx2-1];
				//between region 1 and 2
				else if(abs(X[i])>1e-006 && abs(X[i]-1)>1e-006 && abs(Y[i]-1-dy2)<1e-006)
					T[i]=0.25*(T[i-1]+T[i+1]+T[i-nx1-1]+T[i+nx2+1]);
				//region 2(excluding the boundaries)	
				else
					T[i]=0.25*(T[i-1]+T[i+1]+T[i-nx2-1]+T[i+nx2+1]);						
				}
			}
		}
		
	cout << "iteration complete\n";

	//printing out the result
	ofstream out;
	out.open("proj0002.vtk");
	out << "# vtk DataFile Version 3.1\n";
	out << "temperature distribution\n";
	out << "ASCII\n";
	out << "DATASET UNSTRUCTURED_GRID\n";
	out << "POINTS " << np << " float\n";
	for (int i = 0; i < np; i++)
	{
		out << X[i] << '\t' << Y[i] << '\t' << 0.0 << endl;
	}
	out << "CELLS " << ncell << '\t' << ncell * 5 << endl;
	//region 1
	for(int i=0;i<np1;i++)
	{
		if(abs(X[i]-1.5)<1e-006 || abs(Y[i]-1.0)<1e-006)
			continue;
		else
			out<<4<<'\t'<<i<<'\t'<<i+1<<'\t'<<i+nx1+2<<'\t'<<i+nx1+1<<endl;
	}
	//close to region 2
	for(int i=(nx1+1)*ny1;i<(nx1+1)*ny1+min(nx1,nx2);i++)
	{
			out<<4<<'\t'<<i<<'\t'<<i+1<<'\t'<<i+nx1+2<<'\t'<<i+nx1+1<<endl;
	}
	//region 2
	for(int i=np1;i<np;i++)
	{
		if(abs(X[i]-1.0)<1e-006 || abs(Y[i]-1.5)<1e-006)
			continue;
		else
			out<<4<<'\t'<<i<<'\t'<<i+1<<'\t'<<i+nx2+2<<'\t'<<i+nx2+1<<endl;			
	}
	out << "CELL_TYPES " << np << endl;
	for (int i = 0; i < np; i++)
	{
		out << 9 << endl;
	}
	out << "POINT_DATA " << np << endl;
	out << "SCALARS T float\n";
	out << "LOOKUP_TABLE default\n";
	for (int i = 0; i < np; i++)
	{
		out << T[i] << endl;
	}

	delete[] X;
	delete[] Y;
	delete[] T;
	out.close();
	return 0;
}
