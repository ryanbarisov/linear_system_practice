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
#include <climits>

enum class MatrixFormat
{
	MTX,
	CSR,
	ERROR
};

MatrixFormat GetMatrixFormat(std::string name);

class Matrix
{
public:
	// y <- alpha Ax + beta y
	virtual void Multiply(double alpha, const std::vector<double>& x, double beta, std::vector<double>& y) const = 0;
	virtual int Size() const = 0;
	virtual bool Save(const char * filename, MatrixFormat fmt) const = 0;
	virtual ~Matrix() {}
};


class CSRMatrix : public Matrix
{
private:
	int n;
	bool csr = true;
	std::vector<double> a;
	std::vector<int> ia, ja;
private:
	double eps = 1.0e-12;
	void RemoveZeros(int row,  bool flip_format = false);
	void Add(double alpha, int i, double beta, int j, bool in_place = false, bool flip_format = false);
public:
	CSRMatrix(const std::vector<double>& a, const std::vector<int>& ia, const std::vector<int>& ja);
	int Size() const {return n;}
	double GetEpsilon() const {return eps;}
	void SetEpsilon(double _eps) {eps = _eps;}
	const std::vector<double>& GetA() const {return a;}
	const std::vector<int>& GetIA() const {return ia;}
	const std::vector<int>& GetJA() const {return ja;}
	int GetIA(int row) const {return ia[row];}
	int GetJA(int pos) const {return ja[pos];}
	double GetA(int pos) const {return a[pos];}
	double Get(int i, int j) const;
	void SetA(int pos, double val) {a[pos] = val;}
	void SetIA(int row, int val) {ia[row] = val;}
	void SetJA(int pos, int val) {ja[pos] = val;}
	void FlipStorageFormat(int new_n = -1);
	void Multiply(double alpha, const std::vector<double>& x, double beta, std::vector<double>& y) const;
	void AddRow(double alpha, int i, double beta, int j, bool in_place = false);
	void AddCol(double alpha, int i, double beta, int j, bool in_place = false);
	void RemoveZerosRow(int row);
	void RemoveZerosCol(int col);
	void RemoveZeros();
	bool Save(const char * filename, MatrixFormat fmt) const;
	void PushElement(int row, int col, double elem);
	double Diagonal(int row) const;
	void ExtractDiagonal(std::vector<double>& diag) const;
	double L2Norm(int row) const;
	double L1Norm(int row) const;
};

struct RowAccumulator
{
	int n;
	std::vector<int> jw; // pointer to non-zero elements
	std::vector<double> w; // real values
private:
	std::vector<int> jr; // non-zero indicator
	bool jralloc = false; //
private:
	void split(int p, int pos);
	void split_rec(int p, std::vector<int>& perm, int offset);
	int partition(std::vector<int>& perm, int offset);
	int find_pos(int col) const;
	void remove(int pos);
	bool set(int col, double val, bool add);
	bool prepare_jr();
	bool clear_jr();
public:
	RowAccumulator(int n = 0) : n(n)
	{
		jr.resize(n,-1);
	}
	void Resize(int _n) { if(n != _n) {n = _n; jr.resize(n,-1); jw.clear(); w.clear();}}
	void Clear();
	void SparseAdd(const std::vector<double>& a, const std::vector<int>& ja, double alpha, int jbeg, int jend);
	void SparseAdd(const CSRMatrix* A, int row, double alpha);
	void SetRowFrom(const CSRMatrix* A, int row, int beg = 0, int end = INT_MAX);
	void SetIntervalFrom(int n, const std::vector<int>& rows, const std::vector<double>& vals);
	void SelectLargest(int row, int p);
	void Scale(double alpha);
	void Print() const;
	int Size() const {return w.size();}
	double L2Norm() const;
	double L1Norm() const;
	void Remove(int col);
	double Get(int col) const;
	int GetJW(int pos) const {return jw[pos];}
	double GetW(int pos) const {return w[pos];}
	bool Push(int col, double val) {return set(col, val, false);}
	bool Add(int col, double val) {return set(col, val, true);}
	bool Drop(int col, double norm);
	void Drop(int row, double norm, int p);
	void RemoveZeros(double eps = 1.0e-12);
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

	void push_element(int pos, double val, bool add = false)
	{
		entry e = std::make_pair(pos,val);
		row.push_back(e);
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
	void push_element(int row, int col, double val);

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


int check_nans(const sparse_row& r);
int check_nans(const SparseMatrix* mat);


#endif //MATRIX_H
