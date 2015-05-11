#pragma once

#include "Simplex.h"


class CDVector
{
private:
	bool checkCondition(vect_int& vect, int R);
	void resetValuesBefore(int N, vect_int& vect);
	int  lastNonZero(vect_int& const) const;
	bool checkSolution(matrix_int &m, vect_int& vect, vect_int& centre);
	void considerSolution(vect_int& z, vect_int vect, vect_int& centre);
	inline long calcSolution(vect_int& z, vect_int& vect);
	void BFS(vect_int& z, vect_int& centre, unsigned int R, matrix_int &m);

	vect_int sol_vector;
	long minvalue = 0;
	bool wasFound = false;

public:
	void solveDVector(vect_int& z, matrix_int& conditions, vect_int& initial, unsigned int R);
	void getSolution(vect_int& sol) { sol = sol_vector; }
};

