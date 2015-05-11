#include "Homory.h"


inline double CHomory::fractPart(double fraction)
{
	return (fraction - floor(fraction));
}

int CHomory::getMaxFract()
{
	int res = -1;
	for (unsigned int i = 0; i < m.size()-1; ++i)
	if ((fractPart(m[i][m[i].size() - 1]) > eps) && 
		(fractPart(m[i][m[i].size() - 1]) < 1-eps))
		if ((res == -1) || 
			((res != -1) && (fractPart(m[i][m[i].size() - 1]) > fractPart(m[res][m[res].size() - 1]))))
			res = i;
				
	return res;
}

bool CHomory::initHomory(vect_d &z, matrix &a)
{
	return initSimplex(z, a);
}

bool CHomory::solveHomory()
{
	solveSimplex();
	int res = getMaxFract();

	while (res != -1)
	{
		
		//adding auxiliary inequality
		vect_d aux_ineq = vect_d(m[res].size()+1,0);
		vect_d aux_func(m[res].size() + 1, 0), tmp(m[res].size() + 1, 0);
		double tmp_d = 0;
		
		for (unsigned int i = 0; i < m[res].size() - 1; ++i)
		{
			aux_ineq[i] = fractPart(m[res][i]);
			aux_func[i] = -aux_ineq[i];
		}

		aux_ineq[m[res].size() - 1] = -1;
		aux_func[m[res].size() - 1] = 1;
		aux_ineq[m[res].size()] = fractPart(m[res][m[res].size() - 1]);
		aux_func[m[res].size()] = -fractPart(m[res][m[res].size() - 1]);

		for (int i = 0;  i < m.size(); ++i) // inserting dummy variable for each equation
		{
			tmp_d = m[i][m[i].size() - 1];
			m[i][m[i].size() - 1] = 0;
			m[i].push_back(tmp_d);
		}

		v_free.push_back(v_free.size() + v_base.size());
		v_base.push_back(v_free.size()+v_base.size());
		tmp = m[m.size() - 1];
		m.pop_back();
		m.push_back(aux_ineq);
		m.push_back(tmp);
		m.push_back(aux_func);

		if (!solveSimplex(1)) return false;
		//adding auxiliary inequality

		//ellimination of auxilary variable
		m.pop_back();
		for (int i = 0; i < v_free.size(); ++i)
		if (v_free[i] == v_free.size() + v_base.size() - 1)
		{
			v_free.erase(v_free.begin()+i);
			
			for (int j = 0; j < m.size(); ++j)
				m[j].erase(m[j].begin()+i);

			break;
		}
		if (!solveSimplex()) return false;
		res = getMaxFract();

	} 
	return false;
}

void CHomory::getSolution(vect_d& res)
{
	CSimplex::getSolution(res);
}