#include <method.h>
#include <iomanip>


bool PCG_method::Solve(const std::vector<double>& b, std::vector<double>& x)
{
	int maxiters = params.GetIntegerParameter("maximum_iterations").second;
	double reltol = params.GetRealParameter("relative_tolerance").second;
	double abstol = params.GetRealParameter("absolute_tolerance").second;
	int n = pA->Size();
	std::vector<double> p, Ap(n), r = b, z(n,0.0), z0(n,0.0);
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
		std::copy(z.begin(),z.end(),z0.begin());
		preconditioner->PreconditionedSolve(r,z);
		//Fletcher-Reeves:
		//beta = DotProduct(r,z)/rz;
		//Polak-Ribiere:
		beta = (DotProduct(r,z)-DotProduct(r,z0))/rz;
		Multiply(1,z,beta,p);

		resid = FrobeniusNorm(r);
		iter++;

		// if(iter % 10 == 0)
		{
			std::cout << "iter\t" << std::setw(4) << iter << "\trel_err\t" << std::setw(10) << resid/resid0 << "\t\t|\t" << std::setw(10) << reltol << "\r" << std::flush;
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

bool PBICGStab_method::Solve(const std::vector<double>& b, std::vector<double>& x)
{
	int maxiters = params.GetIntegerParameter("maximum_iterations").second;
	double reltol = params.GetRealParameter("relative_tolerance").second;
	double abstol = params.GetRealParameter("absolute_tolerance").second;
	int n = pA->Size();

	std::vector<double> p(n), r = b, r0, r0_hat, nu0(n,0.0), p0(n, 0.0), y(n, 0.0), h(n), s(n), z(n), t(n);
	pA->Multiply(-1,x,1,r); // r = b-Ax
	r0_hat = r;

	double rho, w, beta;
	double rho0 = 1.0, alpha = 1.0, w0 = 1.0;

	double resid0 = FrobeniusNorm(r), resid;
	if(resid0 < abstol)
	{
		std::cout << "Initial solution satisfies tolerance" << std::endl;
		return true;
	}

	int iter = 0;
	do
	{
		rho = DotProduct(r0_hat, r);
		beta = rho/rho0 * alpha/w0;
		rho0 = rho;
		for(int i = 0; i < n; i++) p[i] = r[i] + beta*(p0[i] - w0*nu0[i]), p0[i] = p[i];
		preconditioner->PreconditionedSolve(p, y);
		pA->Multiply(1.0,y,0.0,nu0); // nu0 <-- Ay
		alpha = rho / DotProduct(r0_hat, nu0);
		h = x;
		Multiply(alpha,y,1.0,h); // h = x + alpha*y
		s = r;
		Multiply(-alpha,nu0,1.0,s); // s = r - alpha*nu0
		preconditioner->PreconditionedSolve(s, z);
		pA->Multiply(1.0,z,0.0,t); 
		w0 = w = DotProduct(t,s)/DotProduct(t,t);
		x = h;
		Multiply(w,z,1.0,x);
		r = s;
		Multiply(-w,t,1.0,r);

		resid = FrobeniusNorm(r);
		iter++;

		// if(iter % 10 == 0)
		{
			std::cout << "iter\t" << std::setw(4) << iter << "\trel_err\t" << std::setw(10) << resid/resid0 << "\t\t|\t" << std::setw(10) << reltol << "\r" << std::flush;
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

// x^new = (I-M^{-1}A) x^old + M^{-1} b
// JACOBI: M = D
void jacobi_precondition(const CSRMatrix* pA, std::vector<double>& x, const std::vector<double>& b)
{
	int n = pA->Size();
	std::vector<double> diag;
	pA->ExtractDiagonal(diag);
	for(int i = 0; i < n; i++)
		x[i] = b[i]/diag[i];
	
}


void gs_precondition(const CSRMatrix* pA, std::vector<double>& x, const std::vector<double>& b)
{
	int n = pA->Size();
	std::vector<double> diag;
	const std::vector<double>& a = pA->GetA();
	const std::vector<int>& ia = pA->GetIA();
	const std::vector<int>& ja = pA->GetJA();
	// (D+L)x = b
	x = b;
	pA->ExtractDiagonal(diag);
	for(int i = 0; i < n; i++)
	{
		for(int j = ia[i]; j < ia[i+1] && ja[j] < i; ++j)
			x[i] -= a[j]*x[ja[j]];
		assert(fabs(diag[i]) > 1.0e-12);
		x[i] /= diag[i];
	}
}

void ssor_precondition(const CSRMatrix* pA, std::vector<double>& x, const std::vector<double>& b)
{
	int n = pA->Size();
	std::vector<double> z = b, diag;
	const std::vector<double>& a = pA->GetA();
	const std::vector<int>& ia = pA->GetIA();
	const std::vector<int>& ja = pA->GetJA();
	pA->ExtractDiagonal(diag);
	// (D+L)z = b
	for(int i = 0; i < n; i++)
	{
		for(int j = ia[i]; j < ia[i+1]; ++j)
			if(ja[j] < i)
				z[i] -= a[j]*z[ja[j]];
		assert(fabs(diag[i]) > 1.0e-12);
		z[i] /= diag[i];
	}
	for(int i = 0; i < n; i++)
		z[i] *= diag[i];
	// (D+U)x = Dz
	x = z;
	for(int i = n-1; i >= 0; i--)
	{
		for(int j = ia[i]; j < ia[i+1]; ++j)
			if(ja[j] > i)
				x[i] -= a[j]*x[ja[j]];
		assert(fabs(diag[i]) > 1.0e-12);
		x[i] /= diag[i];
	}
}

void gs_precondition_backward(const CSRMatrix* pA, std::vector<double>& x, const std::vector<double>& b)
{
	int n = pA->Size();
	std::vector<double> diag;
	const std::vector<double>& a = pA->GetA();
	const std::vector<int>& ia = pA->GetIA();
	const std::vector<int>& ja = pA->GetJA();
	// (D+U)x = b
	x = b;
	pA->ExtractDiagonal(diag);
	for(int i = n-1; i >= 0; i--)
	{
		for(int j = ia[i]; j < ia[i+1]; ++j)
			if(ja[j] > i)
				x[i] -= a[j]*x[ja[j]];
		assert(fabs(diag[i]) > 1.0e-12);
		x[i] /= diag[i];
	}
}

void LU_solve(const CSRMatrix& L, const CSRMatrix& U, const std::vector<double>& b, std::vector<double>& x)
{
	assert(L.Size() == U.Size());
	int n = L.Size();
	assert(b.size() == n && x.size() == n);
	// Lx = b
	{
		const std::vector<double>& a = L.GetA();
		const std::vector<int>& ia = L.GetIA();
		const std::vector<int>& ja = L.GetJA();
		x = b;
		for(int i = 0; i < n; ++i)
			for(int j = ia[i]; j < ia[i+1]; ++j)
				if(i < ja[j])
					x[ja[j]] -= a[j]*x[i];
				//if(ja[j] < i)
					//x[i] -= a[j]*x[ja[j]];
	}
	// Uy = x
	{
		std::vector<double> diag;
		U.ExtractDiagonal(diag);
		const std::vector<double>& a = U.GetA();
		const std::vector<int>& ia = U.GetIA();
		const std::vector<int>& ja = U.GetJA();
		for(int i = n-1; i >= 0; --i)
		{
			for(int j = ia[i]; j < ia[i+1]; ++j)
				if(ja[j] > i)
					x[i] -= a[j]*x[ja[j]];
			assert(fabs(diag[i]) > 1.0e-12);
			x[i] /= diag[i];
		}
	}
}


void LU_in_place_solve(const CSRMatrix* pA, const std::vector<double>& b, std::vector<double>& x)
{
	//assert(A.csr);

	int n = pA->Size();
	const std::vector<int>& ia = pA->GetIA();
	const std::vector<int>& ja = pA->GetJA();
	const std::vector<double>& a = pA->GetA();
	// Lx = b
	x = b;
	for(int i = 0; i < n; ++i)	// loop over rows
	{
		for(int j = ia[i]; j < ia[i+1]; ++j)
		{
			int col = ja[j];
			if(col < i)
				x[i] -= a[j]*x[col];
			else break; // row is sorted by column index in ascending order
		}
	}
	// Uy = x
	std::vector<double> diag;
	pA->ExtractDiagonal(diag);
	for(int i = n-1; i >= 0; --i)	// loop over rows
	{
		for(int j = ia[i+1]-1; j >= ia[i]; --j)
		{
			int col = ja[j];
			if(col > i)
				x[i] -= a[j]*x[col];
			else break; // row is sorted by column index in ascending order
		}
		assert(fabs(diag[i]) > 1.0e-12);
		x[i] /= diag[i];
	}
}
