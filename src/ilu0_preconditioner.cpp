#include <preconditioner.h>
#include <solver_parameters.h>
#include <matrix.h>
#include <method.h>


ILU0_Preconditioner::ILU0_Preconditioner(const CSRMatrix* pA, const SolverParameters& params)
	: Preconditioner(pA, params), LU(nullptr)
	{}

ILU0_Preconditioner::~ILU0_Preconditioner() 
{
	if(LU != nullptr)	delete LU;
}

bool ILU0_Preconditioner::SetupPreconditioner()
{
	LU = new CSRMatrix(*pA);
	construct_inverse(LU);
	return true;
}

bool ILU0_Preconditioner::PreconditionedSolve(const std::vector<double>& rhs, std::vector<double>& x)
{
	LU_in_place_solve(LU,rhs,x);
	return true;
}

void ILU0_Preconditioner::construct_inverse(CSRMatrix* pA)
{
	int n = pA->Size();

	std::vector<double> rA;
	std::vector<int> jA;
	for(int i = 1; i < n; i++)
	{
		jA.clear();
		rA.clear();
		for(int p = pA->GetIA(i), k = pA->GetJA(p); k < i; ++p, k = pA->GetJA(p))
		{
			double aik = pA->GetA(p) / pA->Diagonal(k);
			pA->SetA(p, aik); // A(i,k) <- A(i,k) / A(k,k)
			pA->AddRow(1.0, i, -aik, k, true); // A(i,*) <- A(i,*) - A(i,k)*A(k,*)
			jA.push_back(k);
			rA.push_back(aik);
		}
		// restore A(i,j), j < i
		for(int p = 0; p < jA.size(); ++p)
			pA->PushElement(i, jA[p], rA[p]);
	}
}
