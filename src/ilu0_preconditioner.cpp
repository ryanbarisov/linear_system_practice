#include <preconditioner.h>
#include <solver_parameters.hpp>
#include <matrix.h>
#include <method.h>
#include <list>
#include <matrixwriter.h>


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
	for(int i = 0; i < n; i++)	// loop over rows
	{
		double akk, aik, akj;
		sparse_row& Ai = (*pA)[i];
		for(int j1 = 0; j1 < Ai.row.size(); j1++)
		{
			int k = Ai.row[j1].first;
			if(!(k >= 0 && k < i))	continue;
			akk = diag[k];
			aik = (Ai.row[j1].second /= akk);

			const sparse_row& Ak = (*pA)[k];
			for(int j2 = 0; j2 < Ai.row.size(); j2++)
			{
				int j = Ai.row[j2].first;
				if(!(j >= k+1 && j < n))	continue;
				for(int j3 = 0; j3 < Ak.row.size(); j3++)
				{
					if(Ak.row[j3].first == j)
					{
						akj = Ak.row[j3].second;
						Ai.row[j2].second -= aik*akj;
						if(j == i)
							diag[i] = Ai.row[j2].second;
					}
				}
			}
		}
	}
	
	//MTXMatrixWriter::WriteMatrix(*pA, "test.mtx");
}