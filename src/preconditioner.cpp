#include <preconditioner.h>
#include <method.h>

bool JACOBI_Preconditioner::PreconditionedSolve(const std::vector<double>& rhs, std::vector<double>& x)
{
	jacobi_precondition(pA,x,rhs);
	return true;
}

bool GS_Preconditioner::PreconditionedSolve(const std::vector<double>& rhs, std::vector<double>& x)
{
	gs_precondition(pA,x,rhs); gs_precondition_backward(pA,x,x);
	return true;
}

bool SSOR_Preconditioner::PreconditionedSolve(const std::vector<double>& rhs, std::vector<double>& x)
{
	ssor_precondition(pA,x,rhs);
	return true;
}


PreconditionerType GetPreconditionerType(std::string name)
{
	if (name == "none")
		return PreconditionerType::NONE;
	else if(name == "amg")
		return PreconditionerType::AMG;
	else if(name == "iluc")
		return PreconditionerType::ILUC;
	else if(name == "ilu0")
		return PreconditionerType::ILU0;
	else if(name == "jacobi")
		return PreconditionerType::JACOBI;
	else if(name == "ssor")
		return PreconditionerType::SSOR;
	else if(name == "gs")
		return PreconditionerType::GAUSS_SEIDEL;
	else
	{
		std::cerr << "No such preconditioner: " << name << " don't use preconditioning" << std::endl;
		return PreconditionerType::NONE;
	}
}
