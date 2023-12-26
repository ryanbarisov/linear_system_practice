#include <matrix.h>
#include <matrixwriter.h>


MatrixFormat GetMatrixFormat(std::string filename)
{
	std::string::size_type pos = filename.find_last_of('.');
	if(filename.substr(pos+1) == "mtx")
		return MatrixFormat::MTX;
	else if(filename.substr(pos+1) == "dat")
		return MatrixFormat::CSR;
	else
	{
		std::cerr << "Failed to get matrix format from extension for file: " << filename << std::endl;
		return MatrixFormat::ERROR;
	}
}

void SparseMatrix::add_element(int row, int col, double val)
{
	int N = v.size();
	assert(row >= 0 && row < N);
	entry e = std::make_pair(col, val);
	v[row].add_element(e,true);
}

void SparseMatrix::push_element(int row, int col, double val)
{
	int N = v.size();
	assert(row >= 0 && row < N);
	v[row].push_element(col,val);
}


// y = alpha*A*x + beta*y
void SparseMatrix::Multiply(double alpha, const std::vector<double>& x, double beta, std::vector<double>& y) const
{
	for(int i = 0; i < v.size(); i++)
		y[i] = alpha*v[i].sparse_dot(x) + beta*y[i];
}

bool SparseMatrix::Save(const char * filename, MatrixFormat fmt) const
{
	return WriteMatrix(*this, filename, fmt);
}

void SparseMatrix::Print() const
{
	for(int i = 0; i < v.size(); i++)
	{
		std::cout << "row " << i << std::endl;
		v[i].print();
	}
}

double sparse_dot(const sparse_row& r1, const sparse_row& r2)
{
	double result = 0.0;
	int index1 = 0, index2 = 0, i1 = 0, i2 = 0, column1, column2;
	int size1 = r1.row.size(), size2 = r2.row.size();
	bool process = i1 < size1 && i2 < size2, updated1 = true, updated2 = true;
	while(process)
	{
		if(updated1)	column1 = r1.row[i1].first;
		if(updated2)	column2 = r2.row[i2].first;
		if(column1 < column2)
		{
			i1++;
			updated1 = true;
			updated2 = false;
		}
		else if(column2 < column1)
		{
			i2++;
			updated1 = false;
			updated2 = true;
		}
		else
		{
			result += r1.row[i1].second * r2.row[i2].second;
			i1++;
			i2++;
			updated1 = updated2 = true;
		}
		process = i1 < size1 && i2 < size2;
	}
	return result;
}

double DotProduct(const std::vector<double>& a, const std::vector<double>& b)
{
	assert(a.size() == b.size());
	double prod = 0.0;
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
		y[i] = a*x[i] + b*y[i];
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


int check_nans(const sparse_row& r)
{
	for(int k = 0; k < r.row.size(); k++)
	{
		int col = r.row[k].first;
		double val = r.row[k].second;
		if(std::isnan(val) || std::isinf(val))
		{
			std::cerr << "col " << col << (std::isnan(val) ? " nan " : " inf ") << std::endl;
			return col+1;
		}
	}
	return 0;
}

int check_nans(const SparseMatrix* mat)
{
	int last_found = -1;
	int n = mat->Size();
	for(int i = 0; i < n; i++)
	{
		const sparse_row& r = (*mat)[i];
		if(check_nans(r))
		{
			std::cerr << "(check_nans row " << i << ")." << std::endl;
			last_found = i;
		}
	}
	return last_found;
}

CSRMatrix::CSRMatrix(int _n) : n(n) {}

void CSRMatrix::Multiply(double alpha, const std::vector<double>& x, double beta, std::vector<double>& y) const
{
	if(csc)
		for(int k = 0; k < n; ++k)
		{
			double Ax = 0.0;
			for(int j = ia[k]; j < ia[k+1]; ++j)
				Ax += a[j]*x[ja[k]];
			y[k] = alpha*Ax + beta*y[k];
		}
	else
	{
		std::vector<double> z = y;
		std::fill(y.begin(), y.end(), 0.0);
		for(int k = 0; k < n; ++k)
		{
			for(int j = ia[k]; j < ia[k+1]; ++j)
				y[j] += a[j]*x[k];
			y[k] += beta*z[k];
		}
	}
}



bool CSRMatrix::Save(const char * filename, MatrixFormat fmt) const
{
	return true;
}
