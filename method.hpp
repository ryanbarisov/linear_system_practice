#ifndef METHOD_H
#define METHOD_H

#include "matrix.hpp"
#include <iomanip>

class Method
{
protected:
	const Matrix* ptr_A;
	int maxiters;
	double tol;
	const double reltol = 1e-8;
public:
	virtual bool Setup(const Matrix * A) = 0;
	virtual bool Solve(const std::vector<double>& b, std::vector<double>& x) = 0;
	virtual ~Method() {}
};



class CG_method : public Method
{
public:
	bool Setup(const Matrix * A) 
	{
		ptr_A = A;
		maxiters = A->Size();
		return true;
	}
	bool Solve(const std::vector<double>& b, std::vector<double>& x)
	{
		std::vector<double> p, Ap, r = b;
		ptr_A->Multiply(-1,x,1,r);
		Ap.resize(ptr_A->Size());
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
			ptr_A->Multiply(1,p,0,Ap);
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



class PCG_method : public Method
{
private:
	CSR_matrix * M;	// stores L and U factors together
public:
	bool Setup(const Matrix * A) 
	{
		ptr_A = A;
		M = new CSR_matrix(&*A);
		maxiters = A->Size();
		return true;
	}

	bool Solve(const std::vector<double>& b, std::vector<double>& x)
	{
		std::vector<double> p, Ap, r = b;
		ptr_A->Multiply(-1,x,1,r);
		Ap.resize(ptr_A->Size());
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
			ptr_A->Multiply(1,p,0,Ap);
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


#endif