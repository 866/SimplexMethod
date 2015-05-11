#include <iostream>
#include <time.h>
#include <math.h>
#include "Simplex.h"
#include "DVector.h"
#include "BnB.h"
#include "Simplex.h"
#include "Homory.h"

#define N 10000

using namespace std;

void main()
{
	// solving methods
	CSimplex smplx;
	CHomory hmry;
	CDVector cdv;
	CBnB bnb;

	// performance measuring
	clock_t start;
	double  timeElapsed;

	// Initialize conditions
	matrix_int m(3, vect_int(3, 0));
	vect_int z(2, 0), init(2, 0);
	matrix m_d(2, vect_d(3, 0));
	vect_d z_d(2, 0), sol_d;

	z_d[0] = -350;
	z_d[1] = -150;

	m_d[0][0] = 25;
	m_d[0][1] = 10;
	m_d[0][2] = 100;

	m_d[1][0] = 40;
	m_d[1][1] = 20;
	m_d[1][2] = 190;

	z[0] = -350;
	z[1] = -150;

	m[0][0] = 25;
	m[0][1] = 10;
	m[0][2] = 100;

	m[1][0] = 40;
	m[1][1] = 20;
	m[1][2] = 190;

	init[0] = 0;
	init[1] = 0;
	
	// Solving the task
	cout << N << " iterations";
	start = clock();
	for (int i = 0; i < N; ++i)
	{
		smplx.initSimplex(z_d, m_d);
		smplx.solveSimplex();
		smplx.getSolution(sol_d);
	}
	timeElapsed = double((clock() - start) / (double)CLOCKS_PER_SEC);
	cout << "\nSimplex method. Time elapsed: " << timeElapsed << "\n";
	cout << "Solution: x1 = " << sol_d[0] << "; x2 = " << sol_d[1] << "; f(x1,x2) = " << sol_d[2];
	cout << "\n---------------------------------------------------";

	start = clock();
	for (int i = 0; i < N; ++i)
	{
		bnb.solveBnB(z_d, m_d);
		bnb.getSolution(sol_d);
	}
	timeElapsed = double((clock() - start) / (double)CLOCKS_PER_SEC);
	cout << "\nBnB method. Time elapsed: " << timeElapsed << "\n";
	cout << "Solution: x1 = " << sol_d[0] << "; x2 = " << sol_d[1] << "; f(x1,x2) = " << sol_d[2];
	cout << "\n---------------------------------------------------";

	start = clock();
	for (int i = 0; i < N; ++i)
	{
		cdv.solveDVector(z, m, init, 6);
		cdv.getSolution(init);
	}
	timeElapsed = double((clock() - start) / (double)CLOCKS_PER_SEC);
	cout << "\nDescent Vector Method. Time elapsed: " << timeElapsed << "\n";
	cout << "Solution: x1 = " << init[0] << "; x2 = " << init[1] << "; f(x1,x2) = " << init[0] * z_d[0] + init[1] * z_d[1];
	cout << "\n---------------------------------------------------";

	start = clock();
	for (int i = 0; i < N; ++i)
	{
		hmry.initHomory(z_d, m_d);
		hmry.solveHomory();
		hmry.getSolution(sol_d);
	}
	timeElapsed = double((clock() - start) / (double)CLOCKS_PER_SEC);
	cout << "\nHomory method. Time elapsed: " << timeElapsed << "\n";
	cout << "Solution: x1 = " << sol_d[0] << "; x2 = " << sol_d[1] << "; f(x1,x2) = " << sol_d[2];
	cout << "\n---------------------------------------------------";
	
	system("pause");
	return;
}