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
#include <list>


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


struct MTX_matrix : public Matrix
{
	int M;	// rows
	int N;	// columns
	int L;	// nonzeros

	std::vector<double> a;
	std::vector<int> ia;
	std::vector<int> ja;

public:
	MTX_matrix() : Matrix(MatrixFormat::MTX) {}
	MTX_matrix(const MTX_matrix& other) : Matrix(other), 
		M(other.M), N(other.N), L(other.L), a(other.a), ia(other.ia), ja(other.ja) {}
	MTX_matrix(const char * filename);
	~MTX_matrix() {}

	// y = alpha*A*x + b*y
	void Multiply(double alpha, const std::vector<double>& x, double beta, std::vector<double>& y) const override;
	int Size() const override {return M;}
	bool Save(const char * filename) const override;
};

bool read_rhs_from_mtx(std::vector<double>& rhs, const char* filename);




class CSR_matrix : public Matrix
{
private:
	int n;	// rows
	int nnz;
	double * a;
	int * ia;
	int * ja;
	bool valid = true;
private:
	void swap(CSR_matrix& other);
public:
	CSR_matrix() : Matrix(MatrixFormat::CSR) {}
	CSR_matrix(const char * filename);
	CSR_matrix(const CSR_matrix& other);
	CSR_matrix(const MTX_matrix& mat);
	CSR_matrix(const std::vector<double>& _a, const std::vector<int>& _ia, const std::vector<int>& _ja);

	CSR_matrix& operator=(const CSR_matrix& other);

	~CSR_matrix();

	void Multiply(double alpha, const std::vector<double>& x, double beta, std::vector<double>& y) const override;
	int Size() const override {return n;}
	bool Save(const char * filename) const override;

	// computes x from LUx = b 
	// A stores L and U factors together with the U diagonal inverted
	void LU_solve(const std::vector<double>& b, std::vector<double>& x);
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


class sparse_matrix
{
protected:
	std::vector<sparse_row> v;
public:
	sparse_matrix() : {}
	sparse_matrix(int N) : v(N) {}
	sparse_matrix(const sparse_matrix& other) : v(other.v) {}
	~sparse_matrix() {}

	//void resize(int new_size) {v.resize(new_size);}

	virtual void add_element(int row, int col, double val) = 0;
	virtual void Print() const = 0;

	virtual void Save(const char * filename) const = 0;
	int Size() const {return v.size();}

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

class row_matrix : public sparse_matrix
{
public:
	row_matrix(int N) : sparse_matrix(N) {}
	row_matrix(const row_matrix& other) : sparse_matrix(other) {}
	row_matrix(const MTX_matrix& mat) : sparse_matrix(mat.M)
	{
		for(int i = 0; i < mat.L; i++)
			add_element(mat.ia[i]-1, mat.ja[i]-1, mat.a[i]);
	}

	void add_element(int row, int col, double val) override
	{
		int N = v.size();
		assert(row >= 0 && row < N);
		entry e = std::make_pair(col, val);
		v[row].add_element(e,true);
	}

	// y = alpha*A*x + beta*y
	void Multiply(double alpha, const std::vector<double>& x, double beta, std::vector<double>& y) const override
	{
		for(int i = 0; i < v.size(); i++)
			y[i] = alpha*v[row].sparse_dot(x) + beta*y[i];
	}

	void Print() const override
	{
		for(int i = 0; i < v.size(); i++)
		{
			std::cout << "row " << i << std::endl;
			v[i].print();
		}
	}

	void Save(const char * filename) const override
	{
		std::ofstream ofs(filename);
		ofs << "%% MTX Format matrix" << std::endl;
		for(int i = 0; i < v.size(); i++)
		{
			const sparse_row& r = v[i];
			for(int k = 0; k < r.row.size(); k++)
			{
				int col = r.row[k].first;
				double val = r.row[k].second;
				ofs << i+1 << " " << col+1 << " " << val << std::endl;
			}
		}
		ofs.close();
	}
};


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


#endif