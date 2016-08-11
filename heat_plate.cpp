#include <iostream>
#include <fstream>
using namespace std;

/*------------------------------------------*/
/*source code designed to solve the equation*/
/*of heat transfer using FDM		    */
/*Written by yongbok1102		    */
/*------------------------------------------*/

int main()
{
	//grid
	int nx = 250; int ny = 200;
	double dx = 0.25 / nx; double dy = 0.20 / ny;

	//number of points and cells
	int np = (nx + 1)*(ny + 1);
	int ncell = nx*ny;

	//physical properties
	double k = 122; double h = 50; double Tinf = 30;


	double* X;
	double* Y;
	double* T;
	X = new double[np + 1];
	Y = new double[np + 1];
	T = new double[np + 1];

	int n = 0;

	//generates the grid points
	for (int j = 0; j <= ny; j++)
	{
		for (int i = 0; i <= nx; i++)
		{
			X[n] = i*dx;
			Y[n] = j*dy;
			n++;
		}
	}
	//Initialize the temperature T (solution of the governing equation)
	for (int i = 0; i < np; i++)
	{
		T[i] = 0;
	}

	//Iterative calculation (Gauss-Seidel method)
	for (int itr = 1; itr <= 100000; itr++)
	{
		cout << "iteration " << itr << endl;
		for (int i = 0; i < np; i++)
		{
			//boundary condition
			if (i / (nx + 1) == ny)
			{
				T[i] = (dy*(h / k)*Tinf + T[i - nx - 1]) / (1 + dy*h / k);
			}
			else if (i / (nx + 1) == 0)
			{
				T[i] = T[i + nx + 1];
			}
			else if (i % (nx + 1) == 0)
			{
				T[i] = 50;
			}
			else if (i % (nx + 1) == nx)
			{
				T[i] = 50;
			}
			else
			{
				T[i] = (1. / 4)*(T[i - 1] + T[i + 1] + T[i - nx - 1] + T[i + nx + 1]);
			}
		}
	}
	cout << "iteration complete\n";

	//printing out the result
	ofstream out;
	out.open("proj0001.vtk");
	out << "# vtk DataFile Version 3.1\n";
	out << "temperature distribution\n";
	out << "ASCII\n";
	out << "DATASET UNSTRUCTURED_GRID\n";
	out << "POINTS " << np << " float\n";

	for (int i = 0; i < np; i++)
	{
		out << X[i] << '\t' << Y[i] << '\t' << 0 << endl;
	}
	out << "CELLS " << ncell << '\t' << ncell * 5 << endl;
	for (int i = 0; i < ny*(nx + 1); i++)
	{
		if (i % (nx + 1) == nx) { continue; }
		else {
			out << 4 << '\t' << i << '\t' << i + 1 << '\t' << i + nx + 2 << '\t' << i + nx + 1 << endl;
		}
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
