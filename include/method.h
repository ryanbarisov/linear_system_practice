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
	const CSRMatrix* pA;
public:
	virtual bool Setup(const CSRMatrix * _pA, const SolverParameters& _params)
	{
		pA = _pA;
		params = _params;
		preconditioner = CreatePreconditioner(pA, params);
		return preconditioner ? preconditioner->SetupPreconditioner() : false;
	}
	virtual bool Solve(const std::vector<double>& b, std::vector<double>& x) = 0;
	virtual ~Method() {if(preconditioner) delete preconditioner;}
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


void LU_solve(const CSRMatrix& L, const CSRMatrix& U, const std::vector<double>& b, std::vector<double>& x);
void LU_in_place_solve(const CSRMatrix* pA, const std::vector<double>& b, std::vector<double>& x);

void jacobi_precondition(const CSRMatrix* pA, std::vector<double>& x, const std::vector<double>& b);
void gs_precondition(const CSRMatrix* pA, std::vector<double>& x, const std::vector<double>& b);
void gs_precondition_backward(const CSRMatrix* pA, std::vector<double>& x, const std::vector<double>& b);
void ssor_precondition(const CSRMatrix* pA, std::vector<double>& x, const std::vector<double>& b);

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
