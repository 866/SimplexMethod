#pragma once

#include "Simplex.h"


class CHomory : protected CSimplex
{
private:
	const double eps = 0.000000000001;
	int getMaxFract();
	inline double fractPart(double fraction);

public:
	bool initHomory(vect_d &z, matrix &a);
	bool solveHomory();
	void getSolution(vect_d& res);
};

