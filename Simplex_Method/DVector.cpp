#include "DVector.h"


int CDVector::lastNonZero(vect_int& vect) const
{
	for (int i = vect.size()-1; i >= 0; --i)
	if (vect[i] != 0) return i;
	return 0;
}

bool CDVector::checkSolution(matrix_int &m, vect_int& vect, vect_int& centre) //calc and check conditions
{
	long int sum = 0; 
	for (unsigned int i = 0; i < m.size(); ++i)
	{
		sum = 0;
		for (unsigned int j = 0; j < m[i].size() - 1; ++j)
		{
			if (vect[j] + centre[j]<0) return false;//coords should be non-negative
			sum += m[i][j] * (vect[j] + centre[j]);
		}
		if (sum > m[i][m[i].size() - 1]) return false;
	}
	return true;
}

void CDVector::BFS(vect_int& z, vect_int& centre, unsigned int R, matrix_int &m)
{
	long int sum = 0;
	int i = 0;
	vect_int vect(centre.size(), 0);
	vect[0] = -static_cast<int>(R);
	bool stepBack = true;

	if (checkSolution(m, vect, centre))
		considerSolution(z, vect, centre);

	do
	{
		if ((stepBack) || (sum == 0))
		{
			if (!stepBack)
			{
				if (checkSolution(m, vect, centre))
					considerSolution(z, vect, centre);
			}
			vect[i]++;
			sum -= (vect[i] + vect[i] - 1);
			stepBack = false;
		}
		else
		if (sum < 0)
		{
			sum += (vect[i]*vect[i]);
			vect[i] = 0;
			i--;//step back
			stepBack = true;
		}
		else if (sum > 0)
		{
			if (!stepBack)
			if (vect[i] != 0)
			{
				if (checkSolution(m, vect, centre))
					considerSolution(z, vect, centre);
			}
			else if (lastNonZero(vect) > i)
			{
				if (checkSolution(m, vect, centre))
					considerSolution(z, vect, centre);
			}
			else stepBack = false;

			if (i != vect.size() - 1)
			{
				i++;
				vect[i] = static_cast<int>( -sqrt(sum) );
				sum -= abs(vect[i]*vect[i]);
			}
			else
			{
				vect[i]++;
				sum -= (vect[i] + vect[i] - 1);
			}

		}
		
	} while (vect[0] != R);

	if (checkSolution(m, vect, centre))
		considerSolution(z, vect, centre);

}

void CDVector::considerSolution(vect_int& z, vect_int vect, vect_int& centre)
{
	long int value = 0; // value of function

	for (unsigned int i = 0; i < vect.size(); ++i)
	{
		vect[i] += centre[i]; //displacement
		value += vect[i] * z[i];
	}
	
	if (((wasFound) && (value < minvalue)) || (!wasFound))
	{
		wasFound = true;
		sol_vector = vect;
		minvalue = value;
	}
	
}

void CDVector::solveDVector(vect_int& z, matrix_int& conditions, vect_int& initial, unsigned int R)
{
	sol_vector = initial;
	do
	{
		wasFound = false;
		initial = sol_vector;
		BFS(z, initial, R, conditions);
		if (R>10) R /= 1.1;
	} while (calcSolution(z, initial) - minvalue > 0);
	
	sol_vector = initial;

}

inline long CDVector::calcSolution(vect_int& z, vect_int& vect)
{
	long sum = 0;
	for (int i = 0; i < z.size(); ++i)
		sum += z[i] * vect[i];
	return sum;
}

