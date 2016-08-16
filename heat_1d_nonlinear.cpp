#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
/*---------------------------------------------	*/
/*source code for nonlinear 1D heat conduction 	*/
/*written by yongbok1102			*/
/*---------------------------------------------	*/

//heat conductivity
double k(double T)
{
	static double k0 = 20;
	static double b = 0.00165;
	return k0*(1 + b*T);
}

//residual
double Res(double X[], double Xold[], int N)
{
	double sum = 0;
	for (int i = 0; i<N; i++)
	{
		sum += pow(X[i] - Xold[i], 2);
	}
	return sqrt(sum / N);
}

int main()
{
	//parameters
	double sigma = 5.670*pow(10, -8);
	double L = 0.4;
	double Tinf = 440;

	//mesh generation
	int N = 100; double dx = L / N;
	double* X; X = new double[N + 1];
	for (int i = 0; i <= N; i++)
	{
		X[i] = i*dx;
	}


	double* T; T = new double[N + 1];

	//initialization
	for (int i = 0; i <= N; i++)
	{
		T[i] = 100;
	}
	double* Told; Told = new double[N + 1];

	int itr = 0;

	ofstream out;
	out.open("log.dat");

	//initialization of residual
	double res = 1;

	//Calculated iteratively and printing out the residual
	while (res>1e-013)
	{
		for (int i = 0; i <= N; i++)
		{
			Told[i] = T[i];
		}

		//Gauss-Seidel method for linearized system of equations
		for (int j = 1; j <= 500000; j++)
		{
			for (int i = 0; i <= N; i++)
			{
				if (i == 0)
				{
					T[i] = 460;
				}
				else if (i == N)
				{
					T[i] = (k(Told[i])*(T[i - 1] / dx) + Tinf*sigma*(Told[i] + Tinf)*(pow(Told[i], 2) + pow(Tinf, 2)))
						/ (k(Told[i]) / dx + sigma*(Told[i] + Tinf)*(pow(Told[i], 2) + pow(Tinf, 2)));
				}
				else
				{
					T[i] = (k(0.5*(Told[i + 1] + Told[i]))*T[i + 1] + k(0.5*(Told[i] + Told[i - 1]))*T[i - 1])
						/ (k(0.5*(Told[i + 1] + Told[i])) + k(0.5*(Told[i] + Told[i - 1])));
				}
			}
		}
		itr++;
		res = Res(T, Told, N + 1);
		cout << "iteration: " << itr << " residual: " << res << endl;
		out << itr << '\t' << res << endl;
	}
	out.close();

	//printing out the numerical result
	out.open("res.dat");
	for (int i = 0; i <= N; i++)
	{
		out << X[i] << '\t' << T[i] << endl;
	}
	out.close();

	delete[] X;
	delete[] T;
	delete[] Told;
	return 0;
}
