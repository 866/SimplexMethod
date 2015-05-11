#pragma once

#include <vector>
#include <iostream>


using namespace std;

typedef vector<double> vect_d;
typedef vector<int> vect_int;
typedef vector<vect_int> matrix_int;
typedef vector<vect_d> matrix;

class CSimplex
{
protected:
	vect_int v_base, v_free;
	vect_d coeff;
	matrix m;

	bool isValid(int auxilary = 0) const;
	bool isOptimal(int auxilary = 0) const;
	bool freeIsNonNegative(matrix &a) const;
	void performStep(int auxilary = 0);
	bool auxilaryProblemSolution();

public:
	bool solveSimplex(int auxilary = 0);
	void showTable();
	void getSolution(vect_d &res);
	bool initSimplex(vect_d &z, matrix &a);
};

