#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cassert>
#include <string>
#include <sstream>
#include <iomanip>


double DotProduct(const std::vector<double>& a, const std::vector<double>& b)
{
	assert(a.size() == b.size());
	double prod = 0;
	for(int i = 0; i < a.size(); i++)
		prod += a[i]*b[i];
	return prod;
}

double FrobeniusNorm(const std::vector<double>& v)
{
	return sqrt(DotProduct(v,v));
}

// y = a*x + b*y
void Multiply(double a, const std::vector<double>& x, double b, std::vector<double>& y)
{
	assert(x.size() == y.size());
	for(int i = 0; i < x.size(); i++)
	{
		y[i] = a*x[i] + b*y[i];
	}
}


bool SaveVector(const std::vector<double>& x, const char * filename)
{
	std::ofstream ofs(filename);
	if(!ofs.is_open())
	{
		std::cout << "Failed to open file " << filename << " for writing" << std::endl;
		return false;
	}
	int n = x.size();
	for(int i = 0; i < n; i++)
		ofs << x[i] << std::endl;
	ofs.close();

	return true;
}

enum class MatrixFormat
{
	MTX,
	CSR
};


class Matrix
{
public:
	MatrixFormat fmt;
public:
	Matrix(MatrixFormat _fmt) : fmt(_fmt) {}
	Matrix(const Matrix& other) : fmt(other.fmt) {}
	virtual void Multiply(double alpha, const std::vector<double>& x, double beta, std::vector<double>& y) const = 0;
	virtual int Size() const = 0;
	virtual bool Save(const char * filename) const = 0;
	virtual ~Matrix() {}
};


class MTX_matrix : public Matrix
{
private:
	int M;	// rows
	int N;	// columns
	int L;	// nonzeros

	std::vector<double> a;
	std::vector<int> ia;
	std::vector<int> ja;

	bool valid = true;

	friend class CSR_matrix;
public:
	MTX_matrix() : Matrix(MatrixFormat::MTX) {}
	MTX_matrix(const char * filename) : MTX_matrix()
	{
		std::ifstream ifs(filename);
		if(!ifs.is_open())
		{
			std::cout << "Failed to open file " << filename << std::endl;
			valid = false;
		}
		else
		{
			bool flag = false;
			int i, j, k = 0;
			double aij;
			std::string s;
			while(!ifs.eof())
			{
				std::getline(ifs,s);
				if(s.empty())	continue;
				else if(s[0] == '%')
				{
					continue;	// skip comment 
				}
				else
				{
					std::stringstream ss(s);
					if(!flag)
					{
						ss >> M >> N >> L;
					}
					else
					{
						ss >> i >> j >> aij;
					}
					if(ss.fail())
					{
						std::cout << "Error occurred during reading string: " << s << std::endl;
						valid = false;
						break;
					}
					if(!flag)
					{
						a.resize(L);
						ia.resize(L);
						ja.resize(L);
						flag = true;
					}
					else
					{
						ia[k] = i;
						ja[k] = j;
						a[k] = aij;
						k++;
					}
				}
			}
			assert(k == L);
		}
		ifs.close();
	}
	MTX_matrix(const MTX_matrix& other) : Matrix(other), 
		M(other.M), N(other.N), L(other.L), a(other.a), ia(other.ia), ja(other.ja) {}

	~MTX_matrix() {}

	// y = alpha*A*x + b*y
	void Multiply(double alpha, const std::vector<double>& x, double beta, std::vector<double>& y) const override
	{
		assert(y.size() == M && x.size() == N);
		for(int i = 0; i < M; i++)
		{
			y[i] *= beta;
		}
		for(int k = 0; k < L; k++)
		{
			//int row = ia[k] - 1;
			//int col = ja[k] - 1;
			y[ia[k]-1] += alpha * a[k] * x[ja[k]-1];
		}
	}

	int Size() const override {return M;}

	bool Save(const char * filename) const override
	{
		if(!valid)	return false;
		std::ofstream ofs(filename);
		if(!ofs.is_open())
		{
			std::cerr << "Failed to open file " << filename << " for writing" << std::endl;
			return false;
		}
		ofs << "%%MatrixMarket matrix coordinate real general\n";
		ofs << "% Comment line starts with % sign \n";
		ofs << "% Format: \n";
		ofs << "% \tM N L\n";
		ofs << "% \tI1 J1 A(I1,J1)\n";
		ofs << "% \tI2 J2 A(I2,J2)\n";
		ofs << "% \t ...\n";
		ofs << "% \tIL JL A(IL,JL)\n";
		ofs << "% Indices are 1-based, i.e. A(1,1) is the first element.\n";
		ofs << "%========================================================";

		ofs << M << " " << N << " " << L << std::endl;
		for(int i = 0; i < L; i++)
		{
			ofs << std::setw(7) << ia[i] << "\t" << std::setw(7) << ja[i] << "\t" << std::setw(10) << a[i] << std::endl;
		}


		ofs.close();
		return true;
	}
};


class CSR_matrix : public Matrix
{
private:
	int n;	// rows
	int nnz;
	double * a;
	int * ia;
	int * ja;
	bool valid = true;
public:
	CSR_matrix() : Matrix(MatrixFormat::CSR) {}
	CSR_matrix(const char * filename) : CSR_matrix()
	{
		std::ifstream ifs(filename);
		if(!ifs.is_open())
		{
			std::cout << "Failed to open file " << filename << std::endl;
			valid = false;
		}
		else
		{
			ifs >> n;
			ia = new int[n+1];
			for(int i = 0; i < n+1; i++)
				ifs >> ia[i];
			nnz = ia[n] - ia[0];
			a = new double[nnz];
			ja = new int[nnz];
			for(int i = 0; i < nnz; i++)
			{
				ifs >> ja[i];
			}
			for(int i = 0; i < nnz; i++)
				ifs >> a[i];

			ifs.close();
		}
	}
	CSR_matrix(const CSR_matrix& other) : Matrix(other), n(other.n), nnz(other.nnz)
	{
		if(other.valid)
		{
			a = new double[nnz];
			ia = new int[n+1];
			ja = new int[nnz];
			for(int i = 0; i < nnz; i++)
			{
				a[i] = other.a[i];
				ja[i] = other.ja[i];
			}
			for(int i = 0; i < n+1; i++)
				ia[i] = other.ia[i];
		}
	}
	CSR_matrix(const MTX_matrix& mat) : CSR_matrix()
	{
		if(mat.valid)
		{
			assert(mat.M == mat.N);
			nnz = mat.L;
			n = mat.M;
			
			a = new double[nnz];
			ia = new int[n+1];
			ja = new int[nnz];

			ia[0] = 1;
			for(int i = 0; i < nnz; i++)
			{
				a[i] = mat.a[i];
				ja[i] = mat.ja[i];
				ia[mat.ia[i]]++;
			}
			for(int i = 0; i < n; i++)
				ia[i+1] += ia[i];
		}
	}

	void swap(CSR_matrix& other)
	{
		std::swap(a,other.a);
		std::swap(ia,other.ia);
		std::swap(ja,other.ja);
		std::swap(n,other.n);
		std::swap(nnz,other.nnz);
		std::swap(valid,other.valid);
	}

	CSR_matrix& operator=(const CSR_matrix& other)
	{
		CSR_matrix M(other);
		this->swap(M);
		return *this;
	}
	CSR_matrix(const Matrix* A) : CSR_matrix()
	{
		if(A->fmt == MatrixFormat::CSR)
		{
			const CSR_matrix* ptr_A = dynamic_cast<const CSR_matrix*>(A);
			CSR_matrix M(*ptr_A);
			this->swap(M);
		}
		else if(A->fmt == MatrixFormat::MTX)
		{
			const MTX_matrix* ptr_A = dynamic_cast<const MTX_matrix*>(A);
			CSR_matrix M(*ptr_A);
			this->swap(M);
		}
	}
	CSR_matrix(const std::vector<double>& _a, const std::vector<int>& _ia, const std::vector<int>& _ja) : CSR_matrix()
	{
		n = _ia.size() - 1;
		ia = new int[n+1];
		for(int i = 0; i < n+1; i++)
			ia[i] = _ia[i];
		nnz = _a.size();
		a = new double[nnz];
		ja = new int[nnz];
		for(int i = 0; i < nnz; i++)
		{
			a[i] = _a[i];
			ja[i] = _ja[i];
		}

	}
	~CSR_matrix()
	{
		if(valid)
		{
			delete[] a;
			delete[] ia;
			delete[] ja;
		}
	}

	// y = alpha*A*x + b*y
	void Multiply(double alpha, const std::vector<double>& x, double beta, std::vector<double>& y) const override
	{
		for(int i = 0; i < n; i++)
		{
			y[i] *= beta;
			for(int k = ia[i]; k < ia[i+1]; k++)
			{
				y[i] += alpha * a[k-1] * x[ja[k-1]-1];
			}
		}
	}

	int Size() const override {return n;}


	bool Save(const char * filename) const override
	{
		if(!valid)	return false;
		std::ofstream ofs(filename);
		ofs << n << std::endl;
		for(int i = 0; i < n+1; i++)	
			ofs << ia[i] << "\t";
		ofs << std::endl;
		for(int i = 0; i < nnz; i++)
			ofs << ja[i] << "\t";
		ofs << std::endl;
		for(int i = 0; i < nnz; i++)
			ofs << a[i] << "\t";
		ofs << std::endl;

		ofs.close();
		return true;
	}



	// computes x from LUx = b 
	// A stores L and U factors together with the U diagonal inverted
	void LU_solve(const std::vector<double>& b, std::vector<double>& x)
	{
		assert(b.size() == n && x.size() == n);

		// solve Lx = b
		for(int k = 0; k < n; k++)
		{
			x[k] = b[k];
			// go over nonzero elements of row number k
			for(int i = ia[k]; i < ia[k+1]; i++)
			{
				int col = ja[i-1] - 1;
				// access L factor
				if(col >= k)	break;
				x[k] -= a[i-1] * x[col];
			}
		}

		// solve Uy = x
		for(int k = n-1; k >= 0; k--)
		{
			// go over nonzero elements of row number k
			for(int i = ia[k]; i < ia[k+1]; i++)
			{
				int col = ja[i-1] - 1;
				// access U factor
				if(col <= k)	break;
				x[k] -= a[i-1] * x[col];
			}
		}
		// compute x = x / U(i,i)
		for(int k = 0; k < n; k++)
		{
			for(int i = ia[k]; i < ia[k+1]; i++)
			{
				int col = ja[i-1] - 1;
				if(col == k)	x[k] *= a[i-1];
			}
		}
	}
};


#endif