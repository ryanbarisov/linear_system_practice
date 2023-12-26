#ifndef METHOD_H
#define METHOD_H

#include <matrix.h>
#include <preconditioner.h>
#include <vector>
#include <utility>

enum class MethodType
{
	CG,
	BICGSTAB
};

class Method
{
protected:
	SolverParameters params;
	Preconditioner* preconditioner;
	const SparseMatrix* pA;
public:
	virtual bool Setup(const SparseMatrix * _pA, const SolverParameters& _params)
	{
		pA = _pA;
		params = _params;
		preconditioner = CreatePreconditioner(pA, params);
		return preconditioner ? preconditioner->SetupPreconditioner() : false;
	}
	virtual bool Solve(const std::vector<double>& b, std::vector<double>& x) = 0;
	virtual ~Method() {}
};


class PCG_method : public Method
{
public:
	bool Solve(const std::vector<double>& b, std::vector<double>& x);
};

class PBICGStab_method : public Method
{
public:
	bool Solve(const std::vector<double>& b, std::vector<double>& x);
};


void jacobi_solve(const SparseMatrix* pA, std::vector<double>& x, const std::vector<double>& b);
void gs_solve(const SparseMatrix* pA, std::vector<double>& x, const std::vector<double>& b);

void LU_solve(const SparseMatrix& L, const SparseMatrix& U, const std::vector<double>& b, std::vector<double>& x);
void LU_in_place_solve(SparseMatrix* pA, const std::vector<double>& b, std::vector<double>& x);

void jacobi_precondition(const SparseMatrix* pA, std::vector<double>& x, const std::vector<double>& b);
void gs_precondition(const SparseMatrix* pA, std::vector<double>& x, const std::vector<double>& b);
void symm_gs_precondition(const SparseMatrix* pA, std::vector<double>& x, const std::vector<double>& b);
void gs_precondition_backward(const SparseMatrix* pA, std::vector<double>& x, const std::vector<double>& b);

void sor_precondition(const SparseMatrix* pA, std::vector<double>& x, const std::vector<double>& b);

static Method* CreateMethod(const SolverParameters& params)
{
	Method* method;
	std::string stype = params.GetStringParameter("solver").second;
	if(stype == "cg")
		return new PCG_method();
	else if(stype == "bicgstab")
		return new PBICGStab_method();
	else return NULL;
}

#endif // METHOD_H
