#include <preconditioner.h>
#include <method.h>

bool JACOBI_Preconditioner::PreconditionedSolve(const std::vector<double>& rhs, std::vector<double>& x)
{
	jacobi_precondition(pA,x,rhs);
	return true;
}

bool GS_Preconditioner::PreconditionedSolve(const std::vector<double>& rhs, std::vector<double>& x)
{
	// gauss_seidel_precondition(pA,x,rhs);
	// gauss_seidel_precondition_backward(pA,x,rhs);
	symm_gs_precondition(pA,x,rhs);
	// sor_precondition(pA,x,rhs);
	return true;
}