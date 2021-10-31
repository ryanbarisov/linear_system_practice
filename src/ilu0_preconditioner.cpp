#include <preconditioner.h>
#include <solver_parameters.hpp>
#include <matrix.h>
#include <method.h>
#include <list>


ILU0_Preconditioner::ILU0_Preconditioner(const SparseMatrix* pA, SolverParameters parameters) 
	: Preconditioner(PreconditionerType::ILUC, pA), LU(nullptr) 
	{}

ILU0_Preconditioner::~ILU0_Preconditioner() 
{
	if(LU != nullptr)	delete LU;
}

bool ILU0_Preconditioner::SetupPreconditioner()
{
	int n = pA->Size();
	LU = new SparseMatrix(*pA);
	
	construct_inverse(LU);

	return true;
}

bool ILU0_Preconditioner::PreconditionedSolve(const std::vector<double>& rhs, std::vector<double>& x)
{
	LU_in_place_solve(LU,rhs,x);
	return true;
}

void ILU0_Preconditioner::construct_inverse(SparseMatrix* pA)
{
	int n = pA->Size();

	std::vector<double> diag(n);
	for(int row = 0; row < n; row++)
	{
		const sparse_row& r = (*pA)[row];
		for(int k = 0; k < r.row.size(); k++)
		{
			int col = r.row[k].first;
			if(row == col)
			{
				diag[row] = r.row[k].second;
				break;
			}
		}
	}
	// ilu0
	for(int row = 0; row < n; row++)	// loop over rows
	{
		double akk, aik;
		int cbegin = 0, cend = row;
		sparse_row& r = (*pA)[row];
		for(int j1 = 0; j1 < r.row.size(); j1++)
		{
			int col1 = r.row[j1].first;
			if(!(col1 >= cbegin && col1 > cend))	continue;
			r.row[j1].second /= diag[col1];
			aik = r.row[j1].second;

			const sparse_row& r2 = (*pA)[col1];
			int ccbegin = col1+1, ccend = n;
			for(int j2 = j1+1; j2 < r.row.size(); j2++)
			{
				int col2 = r.row[j2].first;
				if(!(col2 >= ccbegin && col2 < ccend))	continue;
				for(int j3 = 0; j3 < r2.row.size(); j3++)
				{
					if(r2.row[j3].first == col2)
					{
						r.row[j2].second -= aik*r2.row[j3].second;
					}
				}
			}
		}
	}
}