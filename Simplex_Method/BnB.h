#pragma once

#include <math.h>
#include "simplex.h"


class CBnB// Branch and Bound method
{
private:
	struct inequality
	{
		double value = 0;
		unsigned int string = 0;
	};
	typedef vector<vector<inequality>> vec_ineq;
	const double eps = 0.00001;

	int isInteger(vect_d& vect);//returns the number of the element in the vector that is not integer
	int addIneq(vec_ineq& v_i, int num, double value, int string, bool type); // type = true means <=
																			  // type = false means >=
																			  // returns -1 if zero-set
																			  // otherwise - number of string
public:
	void solveBnB(vect_d &z, matrix &a);
	void getSolution(vect_d& v_sol) const { v_sol = min; }

	static vect_d min;
};

