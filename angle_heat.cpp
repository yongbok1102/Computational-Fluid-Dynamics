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

/*----------------------------------------------------------------------------------------------------------*/
/*source code for steady-state heat transfer through the 13cm by 13cm steel angle						   	*/
/*Its thickness is 1.59cm and the boundary conditions are given below									   	*/
/*   																										*/
/*																											*/
/*		       		* *																						*/
/*					* *																						*/
/*		       		* *																						*/
/*					* *																						*/
/*		       		* *																						*/
/*	        600K 	* *																						*/
/*		       		* *																						*/
/*					* *		Exposed to the surrounding air													*/
/*		       		* *		at 300K																			*/
/*					* *		(h = 45W/m2 K)																	*/
/*		       		* *																						*/
/*					* *	* * * * * * * * * * *																*/
/*					* *	* * * * * * * * * * *																*/
/*								600K																		*/
/*Heat conductivity of the steel is given 42.9W/m K															*/
/*----------------------------------------------------------------------------------------------------------*/

int main()
{
	//heat conductivity
	double k = 42.9;
	//convective heat transfer coefficient and temperature of surrounding air
	double h = 45; double Tinf = 300;

	//generating the grid
	int nx1 = 1300; int ny1 = 159;
	int nx2 = 159; int ny2 = 1300 - ny1;
	double dx1 = 0.13 / nx1; double dy1 = 0.0159 / ny1;
	double dx2 = dx1; double dy2 = dy1;
	int np1 = (1 + nx1)*(1 + ny1);
	int np2 = (1 + nx2)*ny2;
	int np = np1 + np2;
	int ncell = nx1*ny1 + nx2*ny2;

	double* X; X = new double[np];
	double* Y; Y = new double[np];
	double* T; T = new double[np];

	int n = 0;
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
			Y[n] = 0.0159 + j*dy2;
			n++;
		}
	}

	//initialization
	for (int i = 0; i < np; i++)
	{
		T[i] = 0;
	}

	//Iterative calculation
	for (int itr = 1; itr <= 20000; itr++)
	{
		cout << "iteration " << itr << endl;


		for (int i = 0; i <np; i++)
		{
			//region 1
			if (i<np1)
			{
				//left side (BC)
				if (i%(nx1+1)==0)
					T[i] = 600;
				//right side (BC)
				else if (i%(nx1+1)==nx1)
					T[i] = (T[i - 1] + (h*dx1*Tinf/k))/(1+(dx1*h/k));
				//lower side (BC)
				else if (i/(nx1+1)==0)
					T[i] = 600;
				//region 1(excluding the boundaries)
				else if (i%(nx1+1)!=0 && i%(nx1+1)!=nx1 && i/(nx1+1)!=0 && i/(nx1+1)!=ny1)
					T[i] = 0.25*(T[i - 1] + T[i + 1] + T[i - nx1 - 1] + T[i + nx1 + 1]);

				//between region 1 and 2
				else
				{	//left side (BC)
					if (i == ny1*(1 + nx1))
						T[i] = 600;
					//upper side(BC)
					else if (i >= ny1*(1 + nx1) + min(nx1, nx2) && i<np1)
						T[i] = (T[i - nx1 - 1] + (h*Tinf*dy1 / k)) / (1 + (dy1*h / k));
					//close to region 2
					else
						T[i] = 0.25*(T[i - 1] + T[i + 1] + T[i - nx1 - 1] + T[i + nx1 + 1]);
				}
			}
			//region 2
			else
			{
				//left side(BC)
				if ((i-np1)%(1+nx2)==0)
					T[i] = 600;
				//right side(BC)
				else if ((i-np1)%(1+nx2)==nx2)
					T[i] = (T[i - 1] + (h*dx1*Tinf / k)) / (1 + (dx1*h / k));
				//upper side(BC)
				else if ((i-np1)/(1+nx2)==ny2-1)
					T[i] = (T[i - nx2 - 1] + (h*Tinf*dy2 / k)) / (1 + (dy2*h / k));
				//between region 1 and 2
				else if ((i-np1)%(1+nx2)!=0 &&  (i-np1)%(1+nx2)!=nx2  && (i-np1)/(1+nx2)==0)
					T[i] = 0.25*(T[i - 1] + T[i + 1] + T[i - nx1 - 1] + T[i + nx2 + 1]);
				//region 2(excluding the boundaries)	
				else
					T[i] = 0.25*(T[i - 1] + T[i + 1] + T[i - nx2 - 1] + T[i + nx2 + 1]);
			}
		}
	}

	cout << "iteration complete\n";

	//printing out the result
	ofstream out;
	out.open("angle_heat.vtk");
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
	for (int i = 0; i<np1; i++)
	{
		if (i%(nx1+1)==nx1 || i/(nx1+1)==ny1)
			continue;
		else
			out << 4 << '\t' << i << '\t' << i + 1 << '\t' << i + nx1 + 2 << '\t' << i + nx1 + 1 << endl;
	}
	//close to region 2
	for (int i = (nx1 + 1)*ny1; i<(nx1 + 1)*ny1 + min(nx1, nx2); i++)
	{
		out << 4 << '\t' << i << '\t' << i + 1 << '\t' << i + nx1 + 2 << '\t' << i + nx1 + 1 << endl;
	}
	//region 2
	for (int i = np1; i<np; i++)
	{
		if ((i-np1)%(1+nx2)==nx2 ||  (i-np1)/(1+nx2) == ny2-1)
			continue;
		else
			out << 4 << '\t' << i << '\t' << i + 1 << '\t' << i + nx2 + 2 << '\t' << i + nx2 + 1 << endl;
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