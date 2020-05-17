#ifndef METHOD_H
#define METHOD_H

#include "matrix.hpp"
#include "preconditioner.hpp"
#include <iomanip>

#include <list>
#include <utility>


class Method
{
protected:
	SolverParameters parameters;
	const MTX_matrix* pA;
	int maxiters;
	double reltol, abstol;
public:
	virtual bool Setup(const MTX_matrix * _pA, PreconditionerType ptype = PreconditionerType::NONE) 
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
	const MTX_matrix* GetMatrix() const {return pA;}

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

	bool Setup(const MTX_matrix * A, PreconditionerType ptype) 
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
	virtual bool Setup(const MTX_matrix* pA, PreconditionerType ptype = PreconditionerType::NONE)
	{
		return Method::Setup(pA);
	}
	virtual bool Solve(const std::vector<double>& b, std::vector<double>& x)
	{
		std::vector<double> p, Ap, r = b;
		pA->Multiply(-1,x,1,r);
		Ap.resize(pA->Size());
		p = r;
		double resid0 = FrobeniusNorm(r), resid;
		if(resid0 < 1e-20)
		{
			std::cout << "Initial solution satisfies tolerance" << std::endl;
			return true;
		}
		double alpha,beta,r2;
		int iter = 0;
		do
		{
			pA->Multiply(1,p,0,Ap);
			r2 = DotProduct(r,r);
			alpha = r2 / DotProduct(Ap,p);
			Multiply(alpha,p,1,x);
			Multiply(-alpha,Ap,1,r);
			beta = DotProduct(r,r)/r2;
			Multiply(1,r,beta,p);

			resid = FrobeniusNorm(r);
			iter++;

			// if(iter % 10 == 0)
			{
				std::cout << "iter\t" << std::setw(4) << iter << "\trel_err\t" << std::setw(10) << resid/resid0 << "\t\t|\t" << std::setw(10) << reltol << "\r";
				std::cout.flush();
				// std::cin.ignore();
			}
			
		} while(iter < maxiters && resid > reltol * resid0);
		if(iter == maxiters && resid > reltol * resid0)
		{
			std::cout << "\nMaximum iterations reached: " << maxiters << " converged to " << resid/resid0 << " relative tolerance " << std::endl;
			return false;
		}
		else
		{
			std::cout << "Iterations: " << iter << " converged to " << resid/resid0 << " relative tolerance " << std::endl;
			return true;
		}
	}
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
	virtual bool Setup(const MTX_matrix* pA, PreconditionerType ptype = PreconditionerType::NONE)
	{
		return PreconditionedMethod::Setup(pA,ptype);
	}

	bool Solve(const std::vector<double>& b, std::vector<double>& x)
	{
		std::vector<double> p, Ap, r = b, z;
		int n = pA->Size();
		Ap.resize(n);
		z.resize(n);
		pA->Multiply(-1,x,1,r);
		preconditioner->PreconditionedSolve(r,z);
		p = z;
		double resid0 = FrobeniusNorm(r), resid;
		if(resid0 < abstol)
		{
			std::cout << "Initial solution satisfies tolerance" << std::endl;
			return true;
		}
		double alpha,beta,rz;
		int iter = 0;
		do
		{
			pA->Multiply(1,p,0,Ap);
			rz = DotProduct(r,z);
			alpha = rz / DotProduct(Ap,p);
			Multiply(alpha,p,1,x);
			Multiply(-alpha,Ap,1,r);
			preconditioner->PreconditionedSolve(r,z);
			beta = DotProduct(r,z)/rz;
			Multiply(1,z,beta,p);

			resid = FrobeniusNorm(r);
			iter++;

			// if(iter % 10 == 0)
			{
				std::cout << "iter\t" << std::setw(4) << iter << "\trel_err\t" << std::setw(10) << resid/resid0 << "\t\t|\t" << std::setw(10) << reltol << "\r";
				std::cout.flush();
				// std::cin.ignore();
			}
			
		} while(iter < maxiters && resid > reltol * resid0);
		if(iter == maxiters && resid > reltol * resid0)
		{
			std::cout << "\nMaximum iterations reached: " << maxiters << " converged to " << resid/resid0 << " relative tolerance " << std::endl;
			return false;
		}
		else
		{
			std::cout << "Iterations: " << iter << " converged to " << resid/resid0 << " relative tolerance " << std::endl;
			return true;
		}
	}
};


#endif
