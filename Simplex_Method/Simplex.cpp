#include "Simplex.h"


bool CSimplex::initSimplex(vect_d& z, matrix& a)
{
	v_free.clear();
	v_base.clear();
	m = a;
	coeff = z;
	m.push_back(z);
	m[m.size()-1].push_back(0);//f(Z) = z0 = 0 at first iteration

	if (freeIsNonNegative(m))
	{
		for (unsigned int i = 0; i < m[0].size()-1; ++i)
			v_free.push_back(i);
		for (unsigned int i = v_free.size(); i < m.size() + v_free.size() - 1; ++i)
			v_base.push_back(i);
		return true;
	}
	else if (auxilaryProblemSolution()) return true;
	return false;
}

bool CSimplex::auxilaryProblemSolution()
{
	vect_int auxilary(0);
	matrix newm(m);
	vect_d tmp(0);

	for (unsigned int i = 0; i < m.size(); ++i)// seaching >= conditions
	{
		if (m[i][m[i].size() - 1] < 0)
		{
			auxilary.push_back(i);
			for (int j = 0; j < m[i].size(); ++j)
				m[i][j] = -m[i][j]; // changing the sign because of >=
		}
		tmp.push_back(m[i][m[i].size() - 1]);
		m[i].resize(m[i].size() - 1);
	}

	
	vect_d aux_func(m[0].size()+auxilary.size(), 0);
	
	for (unsigned int i = 0; i < m[0].size() + auxilary.size(); ++i)
		v_free.push_back(i);

	for (unsigned int i = 0, auxNum = 0, b = 0; i < m.size(); ++i) //making new table with auxiliary variables
	{
		for (int j = 0; j < auxilary.size(); ++j)
			m[i].push_back(0);
			
		if ((auxNum != auxilary.size()) && (i == auxilary[auxNum]))
		{
			m[i][m[i].size() - auxilary.size() + auxNum] = -1;
			if (i != m.size() - 1) v_base.push_back(m.size() + v_free.size() - auxilary.size() + auxNum - 1);
			auxNum++;
			for (int j = 0; j < aux_func.size(); ++j)
				aux_func[j] -= m[i][j];
		}
		else if (i != m.size() - 1)
		{
			v_base.push_back(b+v_free.size());
			b++;
		}
		
	}

	aux_func.push_back(0); // free member
	for (unsigned int i = 0, auxNum = 0; i < m.size(); ++i)
	{
		m[i].push_back(tmp[i]);
		if ((auxNum != auxilary.size()) && (auxilary[auxNum] == i))
		{
			aux_func[aux_func.size()-1] += tmp[i];
			auxNum++;
		}

	}

	tmp = m[m.size() - 1];
	m.push_back(aux_func); 
	if (solveSimplex(1))
	{
		//ellimination of the auxiliary variables
		int auxNum = 0;
		for (unsigned int i = 0; i < v_free.size() && (auxNum != auxilary.size()); ++i)
		if (v_free[i] >= v_base.size() + v_free.size() - auxilary.size())
		{
	
			auxNum++;
			v_free.erase(v_free.begin() + i);
			tmp.erase(tmp.begin() + i);
			for (int j = 0; j < v_base.size() + 1; ++j)
				m[j].erase(m[j].begin()+i);
			i--;
		}

		if (auxNum != auxilary.size()) return false;
		m.pop_back();
		if (solveSimplex()) return true;
			else return false;

	}else return false;
}

bool CSimplex::freeIsNonNegative(matrix &a) const
{
	for (unsigned int i = 0; i < a.size(); ++i)
		if (a[i][a[0].size()-1] < 0) return false;
	
	return true;
}

bool CSimplex::isOptimal(int auxilary) const
{
	bool result = true;

	for (unsigned int i = 0; i < v_free.size(); ++i)
		if (m[v_base.size()+auxilary][i] < 0) result = false;

	return result;
}

bool CSimplex::isValid(int auxilary) const
{
	bool result = false;
	
	for (unsigned int i = 0; ((i < v_free.size()) && (result == false)); ++i)
		if (m[v_base.size()+auxilary][i] < 0)
			for (int j = 0; (j < v_base.size()+auxilary) && (result == false); ++j)
				if (m[j][i] > 0) result = true;

	return result;
}

void CSimplex::performStep(int auxilary)
{
	int pc = 0, pr = -1;//pivot column/raw
	int tmp;

	
	//searching for pivot column
	for (unsigned int i = 0; i < v_free.size(); ++i)
		if (m[v_base.size()+auxilary][i] < m[v_base.size()+auxilary][pc]) pc = i;

	//searching for pivot raw
	double ratio = 0;
	for (unsigned int i = 0; i < v_base.size(); ++i)
	{
		if (m[i][pc] > 0)
		{
			ratio = m[i][v_free.size()] / m[i][pc];
			if (pr == -1) pr = i;
			else if (ratio < (m[pr][v_free.size()] / m[pr][pc])) pr = i;
		}
	}

	//transformation of the table
	for (unsigned int i = 0; i <= v_base.size() + auxilary; ++i)
	for (int j = 0; j <= v_free.size(); ++j)
		if ((i != pr) && (j != pc))
			m[i][j] = (m[pr][pc]*m[i][j]-m[pr][j]*m[i][pc]) / m[pr][pc];
	
	for (unsigned int i = 0; i <= v_free.size(); ++i)
		if (i != pc) m[pr][i] /= m[pr][pc];
	
	for (unsigned int i = 0; i <= v_base.size() + auxilary; ++i)
		if (i != pr) m[i][pc] /= (-m[pr][pc]);

	m[pr][pc] = 1 / m[pr][pc];

	//replacing of elements
	tmp = v_base[pr];
	v_base[pr] = v_free[pc];
	v_free[pc] = tmp;
	
}


bool CSimplex::solveSimplex(int auxilary)
{
	while (isValid(auxilary) && (!isOptimal(auxilary)))
		performStep(auxilary);
	
	return (isValid() || isOptimal());
}

void CSimplex::showTable() // debug feature
{
	for (unsigned int i = 0; i <= v_base.size(); ++i)
	{
		for (int j = 0; j <= v_free.size(); ++j)
			cout << m[i][j]<<" ";
		cout << "\n";
	}

	cout << "v_base: "; 
	for (unsigned int i = 0; i < v_base.size(); ++i)
		cout << v_base[i] << " ";
	cout << "\nv_free: ";
	for (unsigned int i = 0; i < v_free.size(); ++i)
		cout << v_free[i] << " ";

}

void CSimplex::getSolution(vect_d& res)
{
	res.clear();
	res.resize(coeff.size()+1,0);
	for (unsigned int i = 0; i < v_base.size(); ++i)
		if (v_base[i] < v_free.size())
		{
			res[v_base[i]] = m[i][v_free.size()];
			res[v_free.size()] += coeff[v_base[i]] * res[v_base[i]];
		}
}

