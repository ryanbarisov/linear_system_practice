#include <method.h>
#include <iomanip>


bool CG_method::Solve(const std::vector<double>& b, std::vector<double>& x)
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
			std::cout << std::endl;// std::cout.flush();
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


bool PCG_method::Solve(const std::vector<double>& b, std::vector<double>& x)
{
	int n = pA->Size();
	std::vector<double> p, Ap(n), r = b, z(n,0.0);
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
			std::cout << "iter\t" << std::setw(4) << iter << "\trel_err\t" << std::setw(10) << resid/resid0 << "\t\t|\t" << std::setw(10) << reltol << "\n";
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

// SMOOTHING: x^new = G x^old + (I-G) A^{-1} b
// x^new = (I-M^{-1}A) x^old + M^{-1} b
// JACOBI: M = w*D
void jacobi_precondition(const SparseMatrix* pA, std::vector<double>& x, const std::vector<double>& b)
{
	int n = pA->Size();
	std::vector<double> x_old = x;

	x = b;
	pA->Multiply(-1.0,x_old,1.0,x); // x <-- b-Ax

	double omega = 1.0; // Jacobi relaxation parameter, which is used for the optimal smoothing (w = 4/5)
	for(int i = 0; i < n; i++)
	{
		const sparse_row& r = (*pA)[i];
		double aij, aii = 0.0;
		for(int k = 0; k < r.row.size(); k++)
		{
			int j = r.row[k].first;
			aij = r.row[k].second;
			if(j == i)
			{
				x[i] /= (aij * omega);
				break;
			}
		}
	}

	for(int i = 0; i < n; i++)
		x[i] = x_old[i] + x[i];
}

void jacobi_precondition1(const SparseMatrix* pA, std::vector<double>& x, const std::vector<double>& b)
{
	int n = pA->Size();
	x = b;
	for(int i = 0; i < n; i++)
	{
		const sparse_row& r = (*pA)[i];
		double aij, aii = 0.0;
		for(int k = 0; k < r.row.size(); k++)
		{
			int j = r.row[k].first;
			aij = r.row[k].second;
			if(j == i)
			{
				aii = aij; break;
			}
		}
		assert(aii != 0.0);
		x[i] /= aii;
	}
}

// SMOOTHING: x^new = G x^old + (I-G) A^{-1} b
// x^new = (I-M^{-1}A) x^old + M^{-1} b
// GAUSS-SEIDEL: M = D+L
void gs_precondition(const SparseMatrix* pA, std::vector<double>& x, const std::vector<double>& b)
{
	int n = pA->Size();
	std::vector<double> x_old = x;
	x = b;
	pA->Multiply(-1.0,x_old,1.0,x); // x <-- b - Ax
	for(int i = 0; i < n; i++)
	{
		const sparse_row& r = (*pA)[i];
		double aij, aii = 0.0;
		for(int k = 0; k < r.row.size(); k++)
		{
			int j = r.row[k].first;
			aij = r.row[k].second;
			if(j == i)
				aii = aij;
			else if (j < i)
				x[i] -= aij * x[j];
		}
		assert(aii != 0.0);
		x[i] /= aii;
	}
	for(int i = 0; i < n; i++)
		x[i] = x_old[i] + x[i];
}

void gs_precondition1(const SparseMatrix* pA, std::vector<double>& x, const std::vector<double>& b)
{
	double w = 1.0;

	int n = pA->Size();
	x = b;
	for(int i = 0; i < n; i++)
	{
		const sparse_row& r = (*pA)[i];
		double aij, aii = 0.0;
		for(int k = 0; k < r.row.size(); k++)
		{
			int j = r.row[k].first;
			aij = r.row[k].second;
			if(j == i)
				aii = aij;
			else if (j < i)
				x[i] -= aij * x[j];
		}
		assert(aii != 0.0);
		x[i] /= aii;
	}
}

void symm_gs_precondition(const SparseMatrix* pA, std::vector<double>& x, const std::vector<double>& b)
{
	double w = 1.0;

	int n = pA->Size();
	std::vector<double> z;
	z = b;
	for(int i = 0; i < n; i++)
	{
		const sparse_row& r = (*pA)[i];
		double aij, aii = 0.0;
		for(int k = 0; k < r.row.size(); k++)
		{
			int j = r.row[k].first;
			aij = r.row[k].second;
			if(j == i)
			{
				aii = aij;
				z[i] *= aii;
				break;
			}
		}
	}
	// (D+L)z = Db
	for(int i = 0; i < n; i++)
	{
		const sparse_row& r = (*pA)[i];
		double aij, aii = 0.0;
		for(int k = 0; k < r.row.size(); k++)
		{
			int j = r.row[k].first;
			aij = r.row[k].second;
			if(j == i)
				aii = aij;
			else if(j < i)
				z[i] -= aij*z[j];
		}
		assert(aii != 0.0);
		z[i] /= aii;
	}
	// (D+U)x = z
	x = z;
	for(int i = n-1; i >= 0; i--)
	{
		const sparse_row& r = (*pA)[i];
		double aij, aii = 0.0;
		for(int k = 0; k < r.row.size(); k++)
		{
			int j = r.row[k].first;
			aij = r.row[k].second;
			if(j == i)
				aii = aij;
			else if(j > i)
				x[i] -= aij*x[j];
		}
		assert(aii != 0.0);
		x[i] /= aii;
	}
}

void gs_precondition_backward(const SparseMatrix* pA, std::vector<double>& x, const std::vector<double>& b)
{
	double w = 1.0;

	int n = pA->Size();
	x = b;
	for(int i = n-1; i >= 0; i--)
	{
		const sparse_row& r = (*pA)[i];
		double aij, aii = 0.0;
		for(int k = 0; k < r.row.size(); k++)
		{
			int j = r.row[k].first;
			aij = r.row[k].second;
			if(j == i)
				aii = aij;
			else if (j > i)
				x[i] -= aij * x[j];
		}
		assert(aii != 0.0);
		x[i] /= aii;
	}
}

void sor_precondition(const SparseMatrix* pA, std::vector<double>& x, const std::vector<double>& b)
{
	int n = pA->Size();
	std::vector<double> x_old(n, 0.0);
	gs_precondition(pA,x,b);
	
	double w = 1.0;
	if(w != 1.0)
		for(int i = 0; i < n; i++)
		{
			x[i] = (1.0-w)*x_old[i] + w*x[i];
		}
}




void jacobi_solve(const SparseMatrix* pA, std::vector<double>& x, const std::vector<double>& b)
{
	int n = pA->Size();
	std::vector<double> x_old = x;
	
	for(int i = 0; i < n; i++)
	{
		const sparse_row& r = (*pA)[i];
		double aij, aii = 0.0;
		x[i] = b[i];
		for(int k = 0; k < r.row.size(); k++)
		{
			int j = r.row[k].first;
			aij = r.row[k].second;
			if(j == i)
				aii = aij;
			else
				x[i] -= aij*x_old[j];
		}
		assert(aii != 0.0);
		x[i] /= aii;
	}
}

void gs_solve(const SparseMatrix* pA, std::vector<double>& x, const std::vector<double>& b)
{
	double w = 1.0;

	int n = pA->Size();
	std::vector<double> x_old = x;
	for(int i = 0; i < n; i++)
	{
		const sparse_row& r = (*pA)[i];
		double aij, aii = 0.0;
		x[i] = b[i];
		for(int k = 0; k < r.row.size(); k++)
		{
			int j = r.row[k].first;
			aij = r.row[k].second;
			if(j == i)
				aii = aij;
			else if(j > i)
				x[i] -= aij * x_old[j];
			else // if (j < i)
				x[i] -= aij * x[j];
		}
		assert(aii != 0.0);
		x[i] /= aii;
	}
	if(w != 1.0)
		for(int i = 0; i < n; i++)
		{
			x[i] = (1.0-w)*x_old[i] + w*x[i];
		}
}


void LU_solve(const SparseMatrix& L, const SparseMatrix& U, const std::vector<double>& b, std::vector<double>& x)
{
	assert(L.Size() == U.Size());
	int n = L.Size();
	assert(b.size() == n && x.size() == n);

	// solve Lx = b
	for(int k = 0; k < n; k++)
	{
		x[k] = b[k];
		const sparse_row& r = L[k];
		for(sparse_type::const_iterator it = r.row.begin(); it != r.row.end(); it++)
		{
			if(it->first < k)
				x[k] -= it->second * x[it->first];
		}
	}

	// solve Uy = x
	for(int k = n-1; k >= 0; k--)
	{
		const sparse_row& r = U[k];
		for(sparse_type::const_iterator it = r.row.begin(); it != r.row.end(); it++)
		{
			if(it->first > k)
				x[k] -= it->second * x[it->first];
		}
		x[k] /= r.get_element(k);
	}
}


void LU_in_place_solve(SparseMatrix* pA, const std::vector<double>& b, std::vector<double>& x)	
{
	int n = pA->Size();

	// solve Lx = b
	x = b;
	for(int i = 0; i < n; i++)	// loop over rows
	{
		const sparse_row& r = (*pA)[i];
		for(int kk = 0; kk < r.row.size(); kk++)
		{
			int j = r.row[kk].first;
			double aij = r.row[kk].second;
			if(j < i)
				x[i] -= aij*x[j];
			else break; // row is sorted by column index in ascending order
		}
	}
	// solve Uy = x
	for(int i = n-1; i >= 0; i--)	// loop over rows
	{
		double aii = 0.0;
		const sparse_row& r = (*pA)[i];
		for(int kk = r.row.size()-1; kk >= 0; kk--)
		{
			int j = r.row[kk].first;
			double aij = r.row[kk].second;
			if(j == i)
				aii = aij;
			else if(j > i)
				x[i] -= aij*x[j];
			else break; // row is sorted by column index in ascending order
		}
		assert(aii != 0.0);
		x[i] /= aii;
	}
}