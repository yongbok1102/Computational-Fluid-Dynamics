#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
/*------------------------------------------*/
/*source code for the potential flow		    */
/*solution method for PDE: FDM				      */
/*iteration method: Gauss-Seidel method		  */
/*Written by L.Y.B							            */
/*------------------------------------------*/

int main()
{
	//grid
	int nx = 100; int ny = 50;
	double dx = 10. / nx; double dy = 5. / ny;
	
	//number of points and cells
	int np = (nx + 1)*(ny + 1);
	int ncell = nx*ny;

	double* X;
	double* Y;
	double* PSI;
	double* U; double* V;
	X = new double[np];
	Y = new double[np];
	PSI = new double[np];	
	U = new double[np];
	V = new double[np];
	
	int n = 0;

	//numbering the grid points
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
		PSI[i] = 0;
	}
	
	//Iterative calculation (Gauss-Seidel method)
	for (int itr = 1; itr <= 10000; itr++)
	{
		cout << "iteration " << itr << endl;
		for (int i = 0; i < np; i++)
		{
			//lower side
			if(abs(Y[i])<1e-006)
				PSI[i]=0;
			//upper side
			else if(abs(Y[i]-dy*ny)<1e-006)
				PSI[i]=5;
			//left side
			else if(abs(X[i])<1e-006 && abs(Y[i])>1e-006 && abs(Y[i]-dy*ny)>1e-006)
				PSI[i]=Y[i];
			//right side
			else if (abs(X[i]-dx*nx)<1e-006 && abs(Y[i])>1e-006 && abs(Y[i]-dy*ny)>1e-006)
			{
				if(i/(nx+1)<=15)
					PSI[i]=0;
				else if(i/(nx+1)>=35)
					PSI[i]=5;
				else
					PSI[i]=PSI[i-1];
			}
			else
				PSI[i]=0.25*(PSI[i-1]+PSI[i+1]+PSI[i-nx-1]+PSI[i+nx+1]);
		}
	}

	//calculating the velocity field
		//X-comp
		for(int i=0;i<np;i++)
		{
			if(abs(Y[i])<1e-006)
				U[i]=(PSI[i+nx+1]-PSI[i])/dy;
			else if(abs(Y[i]-ny*dy)<1e-006)
				U[i]=(PSI[i]-PSI[i-nx-1])/dy;
			else
				U[i]=0.5*(PSI[i+nx+1]-PSI[i-nx-1])/dy;
		}	
		//Y-comp
		for(int i=0;i<np;i++)
		{
			if(abs(X[i])<1e-006)
				V[i]=-(PSI[i+1]-PSI[i])/dx;
			else if(abs(X[i]-nx*dx)<1e-006)
				V[i]=-(PSI[i]-PSI[i-1])/dx;
			else
				V[i]=-0.5*(PSI[i+1]-PSI[i-1])/dx;
		}
	cout << "iteration complete\n";

	//printing out the result
	ofstream out;
	out.open("proj0003.vtk");
	out << "# vtk DataFile Version 3.1\n";
	out << "potential flow field\n";
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
	out << "SCALARS PSI float\n";
	out << "LOOKUP_TABLE default\n";
	for (int i = 0; i < np; i++)
	{
		out << PSI[i] << endl;
	}
	out<<"VECTORS U float\n";
	for(int i=0;i<np;i++)
	{
		out<<U[i]<<'\t'<<V[i]<<'\t'<<0.0<<endl;
	}

	delete[] X;
	delete[] Y;
	delete[] PSI;
	delete[] U;
	delete[] V;
	out.close();
	return 0;	
	
}
