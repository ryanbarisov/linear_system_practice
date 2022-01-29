#include <preconditioner.h>
#include <solver_parameters.hpp>
#include <matrix.h>
#include <method.h>
#include <algorithm>
#include <set>

AMG_Preconditioner::AMG_Preconditioner(const SparseMatrix* _pA, SolverParameters parameters)
	: Preconditioner(PreconditionerType::AMG, pA)
{
	
}