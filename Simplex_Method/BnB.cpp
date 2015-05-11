#include "BnB.h"


vect_d CBnB::min;

int CBnB::isInteger(vect_d &vect)
{
	for (unsigned int i = 0; i < vect.size(); ++i)
		if (fabs(round(vect[i]) - vect[i]) > eps) return i;
	
	return -1;
}

void CBnB::solveBnB(vect_d &z, matrix &a)
{
	vect_d solution;
	CSimplex simplx_meth;

	if (simplx_meth.initSimplex(z, a))
	if (simplx_meth.solveSimplex())
	{
		simplx_meth.getSolution(solution);
		unsigned int isInt = isInteger(solution);

		if (isInt == -1)
		{
			if (((min.size() == 0) && (isInt == -1))
				|| ((min.size() != 0) && (min[min.size() - 1] >= solution[solution.size() - 1])))
				min = solution;
		}
		else if (isInt != (min.size() - 1))
		{
			vect_d new_free(a[0].size(), 0);
			matrix newa(a);
			
			///xb + x(n+1) <= solution
				new_free[isInt] = 1;
				new_free[new_free.size() - 1] = trunc(solution[isInt]);
				newa.push_back(new_free);
				solveBnB(z, newa);//lower bound
	
			
				///xb + x(n+1) >= solution
				if (trunc(solution[isInt]) != 0)
				{
					newa[newa.size() - 1][isInt] = -1;
					newa[newa.size() - 1][new_free.size() - 1] = -trunc(solution[isInt] + 1);
					solveBnB(z, newa);//upper bound
				}

		}
	}
		
}