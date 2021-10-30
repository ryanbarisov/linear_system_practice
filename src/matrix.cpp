#include <matrix.h>
#include <matrixwriter.h>



void SparseMatrix::add_element(int row, int col, double val)
{
	int N = v.size();
	assert(row >= 0 && row < N);
	entry e = std::make_pair(col, val);
	v[row].add_element(e,true);
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