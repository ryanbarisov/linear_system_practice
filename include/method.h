#ifndef METHOD_H
#define METHOD_H

#include <matrix.h>
#include <preconditioner.h>
#include <vector>
#include <utility>


class Method
{
protected:
	SolverParameters parameters;
	const SparseMatrix* pA;
	int maxiters;
	double reltol, abstol;
public:
	virtual bool Setup(const SparseMatrix * _pA, PreconditionerType ptype = PreconditionerType::NONE) 
	{
		pA = _pA;
		maxiters = GetIntegerParameter("maximum_iterations");
		reltol = GetRealParameter("relative_tolerance");
		abstol = GetRealParameter("absolute_tolerance");
		return true;
	}
	virtual bool Solve(const std::vector<double>& b, std::vector<double>& x) = 0;
	virtual ~Method() {}

	Method()
	{
		// common default parameters goes here
		parameters.SetIntegerParameter("maximum_iterations", 2500);
		parameters.SetRealParameter("relative_tolerance", 1.0e-8);
		parameters.SetRealParameter("absolute_tolerance", 1.0e-12);
	}
	const SparseMatrix* GetMatrix() const {return pA;}

	void SetRealParameter(std::string name, double value)	{parameters.SetRealParameter(name, value);}
	void SetIntegerParameter(std::string name, int value)	{parameters.SetIntegerParameter(name, value);}
	double GetRealParameter(std::string name) const
	{
		std::pair<bool,double> find_pair = parameters.GetRealParameter(name);
		if(!find_pair.first)
			std::cerr << "Parameter " << name << " was not set." << std::endl;
		return find_pair.second;
	}
	int GetIntegerParameter(std::string name) const
	{
		std::pair<bool,int> find_pair = parameters.GetIntegerParameter(name);
		if(!find_pair.first)
			std::cerr << "Parameter " << name << " was not set." << std::endl;
		return find_pair.second;
	}
};

class PreconditionedMethod : public Method
{
protected:
	Preconditioner * preconditioner;
public:
	PreconditionedMethod() : preconditioner(NULL)
	{
		// default preconditioner parameters goes here
		parameters.SetRealParameter("drop_tolerance", 1.0e-4);
		parameters.SetIntegerParameter("level_of_fill", 80);
	}
	PreconditionedMethod(const PreconditionedMethod& other)
	{
		preconditioner = CreatePreconditioner(other.preconditioner->GetType(), other.GetMatrix(), other.parameters);
	}
	PreconditionedMethod& operator=(const PreconditionedMethod& other)
	{
		if(preconditioner != NULL)	
		{
			delete preconditioner;
			preconditioner = NULL;
		}
		preconditioner = CreatePreconditioner(other.preconditioner->GetType(), other.GetMatrix(), other.parameters);
		return *this;
	}
	virtual ~PreconditionedMethod()
	{
		if(preconditioner != NULL)	
		{
			delete preconditioner;
			preconditioner = NULL;
		}
	}

	bool Setup(const SparseMatrix * A, PreconditionerType ptype) 
	{
		Method::Setup(A);
		if(preconditioner == NULL)
			preconditioner = CreatePreconditioner(ptype, A, parameters);
		if(preconditioner == NULL)	return false;
		return preconditioner->SetupPreconditioner();
	}
};




class CG_method : public Method
{
public:
	virtual bool Setup(const SparseMatrix* pA, PreconditionerType ptype = PreconditionerType::NONE)
	{
		return Method::Setup(pA);
	}
	bool Solve(const std::vector<double>& b, std::vector<double>& x);
};



class PCG_method : public PreconditionedMethod
{
public:
	PCG_method() : PreconditionedMethod() {}
	PCG_method(const PCG_method& other) : PreconditionedMethod(other) {}
	PCG_method& operator=(const PCG_method& other) 
	{
		PreconditionedMethod::operator=(other);
		return *this;
	}
	virtual ~PCG_method()	{}
	virtual bool Setup(const SparseMatrix* pA, PreconditionerType ptype = PreconditionerType::NONE)
	{
		return PreconditionedMethod::Setup(pA,ptype);
	}

	bool Solve(const std::vector<double>& b, std::vector<double>& x);
};


void jacobi(const SparseMatrix* pA, std::vector<double>& x, const std::vector<double>& b);
void gauss_seidel(const SparseMatrix* pA, std::vector<double>& x, const std::vector<double>& b);
void LU_solve(const SparseMatrix& L, const SparseMatrix& U, const std::vector<double>& b, std::vector<double>& x);


#endif // METHOD_H
