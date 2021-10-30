#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>


enum class MatrixFormat
{
	MTX,
	CSR,
	ERROR
};

class Matrix
{
public:
	virtual void Multiply(double alpha, const std::vector<double>& x, double beta, std::vector<double>& y) const = 0;
	virtual int Size() const = 0;
	virtual bool Save(const char * filename, MatrixFormat fmt) const = 0;
	virtual ~Matrix() {}
};

typedef std::pair<int,double> entry;
typedef std::vector<entry> sparse_type;

struct sparse_row
{
	sparse_type row;

	sparse_row() {}

	// assign part of row in interval [begin, end)
	sparse_row& assign(const sparse_row& other, int begin, int end)
	{
		row.clear();
		if(begin < end)
		{
			for(sparse_type::const_iterator it = other.row.begin(); it != other.row.end(); ++it)
				if(it->first >= begin && it->first < end)
					row.push_back(*it);
		}
		return *this;
	}

	void add_element(const entry& e, bool add = false)
	{
		if(row.empty())	
		{
			row.push_back(e);
			return;
		}
		sparse_type::iterator pos = row.end();
		bool equal = false;
		int r = e.first;
		for(sparse_type::iterator it = row.begin(); it != row.end(); ++it)
		{
			if(r <= it->first)
			{
				pos = it;
				if(r == it->first)	equal = true;
				break;
			}
		}
		if(equal)
		{
			if(add)	pos->second += e.second;
			else pos->second = e.second;
		}
		else if(pos == row.end())
			row.push_back(e);
		else 
			row.insert(pos, e);
	}

	void add_element(int pos, double val, bool add = false)
	{
		entry e = std::make_pair(pos,val);
		add_element(e, add);
	}

	sparse_row& operator+=(const sparse_row& other)
	{
		for(sparse_type::const_iterator it = other.row.begin(); it != other.row.end(); ++it)
			add_element(*it, true);
		return *this;
	}

	sparse_row operator+(const sparse_row& other)
	{
		sparse_row r = *this;
		r += other;
		return r;
	}

	// update row as row + a * other.row in interval [begin,end)
	sparse_row& plus(double a, const sparse_row& other, int begin, int end)
	{
		if(begin < end)
		{
			entry e;
			for(sparse_type::const_iterator it = other.row.begin(); it != other.row.end(); ++it)
			{
				if(it->first >= begin && it->first < end)
				{
					e.first = it->first;
					e.second = it->second * a;
					add_element(e, true);
				}
			}
		}
		return *this;
	}

	sparse_row& operator*=(double alpha)
	{
		for(sparse_type::iterator it = row.begin(); it != row.end(); ++it)
			it->second *= alpha;
		return *this;
	}

	sparse_row operator*(double alpha) const
	{
		sparse_row r = *this;
		r *= alpha;
		return r;
	}

	void remove_zeros() 
	{
		for(int i = 0; i < row.size();)
		{
			if(fabs(row[i].second) < 1e-10)
				row.erase(row.begin()+i);
			else
				i++;
		}
	}

	double get_element(int index) const
	{
		for(sparse_type::const_iterator it = row.begin(); it != row.end(); ++it)
			if(it->first == index)	return it->second;
		return 0.0;
	}

	void print() const
	{
		for(sparse_type::const_iterator it = row.begin(); it != row.end(); ++it)
			std::cout << it->first << ": " << it->second << std::endl;
	}

	double sparse_dot(const std::vector<double>& v) const
	{
		double result = 0.0;
		int nv = v.size();
		for(int k = 0; k < row.size(); k++)
		{
			int col = row[k].first;
			assert(col < nv);
			result += row[k].second * v[col];
		}
		return result;
	}
};


class SparseMatrix : public Matrix
{
protected:
	std::vector<sparse_row> v;
public:
	SparseMatrix(int N) : v(N) {}
	SparseMatrix(const SparseMatrix& other) : v(other.v) {}
	~SparseMatrix() {}

	void add_element(int row, int col, double val);

	// y = alpha*A*x + beta*y
	void Multiply(double alpha, const std::vector<double>& x, double beta, std::vector<double>& y) const override;
	int Size() const override {return v.size();}
	bool Save(const char * filename, MatrixFormat fmt) const override;
	void Print() const;


	void set_vector(int pos, const sparse_row& vec, int begin, int end)
	{
		assert(pos >= 0 && pos < v.size());
		v[pos].assign(vec,begin,end);
	}
	void add_vector(int pos, double coef, const sparse_row& vec, int begin, int end)
	{
		assert(pos >= 0 && pos < v.size());
		v[pos].plus(coef, vec, begin, end);
	}

	const sparse_row& operator[](int pos) const 
	{
		assert(pos >= 0 && pos < v.size());
		return v[pos];
	}

	sparse_row& operator[](int pos) 
	{
		assert(pos >= 0 && pos < v.size());
		return v[pos];
	}
};


double sparse_dot(const sparse_row& r1, const sparse_row& r2);

double DotProduct(const std::vector<double>& a, const std::vector<double>& b);
double FrobeniusNorm(const std::vector<double>& v);
void Multiply(double a, const std::vector<double>& x, double b, std::vector<double>& y);
bool SaveVector(const std::vector<double>& x, const char * filename);


#endif //MATRIX_H