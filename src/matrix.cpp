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

CSRMatrix::CSRMatrix(const std::vector<double>& a, const std::vector<int>& ia, const std::vector<int>& ja) : a(a), ia(ia), ja(ja)
{
	n = ia.size()-1;
}

void CSRMatrix::FlipStorageFormat(int new_n)
{

	if(!a.empty())
	{
		if(new_n == -1) new_n = n;
		std::vector<int> new_ia(new_n+1, 0), new_ja(ja.size());
		std::vector<double> new_a(a.size());
		for(int k = 0; k < n; ++k)
			for(int j = ia[k]; j < ia[k+1]; ++j)
				new_ia[ja[j]+1]++;
		for(int j = 1; j < new_n+1; ++j)
			new_ia[j] += new_ia[j-1];
		std::vector<int> pos = new_ia;
		for(int k = 0; k < n; ++k)
			for(int j = ia[k]; j < ia[k+1]; ++j)
			{
				new_ja[pos[ja[j]]] = k;
				new_a[pos[ja[j]]] = a[j];
				++pos[ja[j]];
			}
		a = new_a;
		ia = new_ia;
		ja = new_ja;
		n = new_n;
	}
	csr = !csr;
}

void CSRMatrix::Multiply(double alpha, const std::vector<double>& x, double beta, std::vector<double>& y) const
{
	if(csr)
		for(int k = 0; k < n; ++k)
		{
			double Ax = 0.0;
			for(int j = ia[k]; j < ia[k+1]; ++j)
				Ax += a[j]*x[ja[j]];
			y[k] = alpha*Ax + beta*y[k];
		}
	else
	{
		std::vector<double> z = y;
		std::fill(y.begin(), y.end(), 0.0);
		for(int k = 0; k < n; ++k)
		{
			for(int j = ia[k]; j < ia[k+1]; ++j)
				y[ja[j]] += a[j]*x[k];
			y[k] += beta*z[k];
		}
	}
}

// A(i,*) <- alpha * A(i,*) + beta * A(j,*)
void CSRMatrix::Add(double alpha, int rowi, double beta, int rowj, bool in_place, bool flip_format)
{
	assert(rowi < n && rowj < n);
	bool format = !flip_format ? csr : !csr;
	if(format)
	{
		int i0 = ia[rowi];
		if(fabs(1.0-alpha) > 1.0e-20)
			for(int j = ia[rowi]; j < ia[rowi+1]; ++j)
				a[j] *= alpha;
		if(beta)
			for(int j = ia[rowj]; j < ia[rowj+1]; ++j)
			{
				double addend = beta*a[j];
				int icol = ja[i0], i1 = ia[rowi+1], jcol = ja[j];
				while(jcol > icol && i0 < i1) icol = ja[++i0];
				if(!in_place && (jcol < icol || i0 == i1))
				{
					a.insert(a.begin()+i0, addend);
					ja.insert(ja.begin()+i0, jcol);
					if(j > i0) ++j;
					for(int k = rowi+1; k < n+1; ++k)
						ia[k]++;
				}
				else if(jcol == icol)
					a[i0] += addend;
			}
	}
	else
	{
		std::vector<int> rows, cols;
		std::vector<double> elems;
		rows.reserve(n);
		elems.reserve(n);
		for(int k = 0; k < n; ++k)
			for(int j = ia[k]; j < ia[k+1]; ++j)
				if(ja[j] == rowi)
					a[j] *= alpha;
				else if(beta && ja[j] == rowj)
				{
					rows.push_back(ja[j]);
					cols.push_back(k);
					elems.push_back(a[j]);
				}
		if(!cols.empty())
		{
			int j0 = 0, j1 = cols.size();
			for(int k = 0; k < n && j0 < j1; ++k) if(k == cols[j0])
			{
				bool found = false;
				int pos = ia[k];
				for(int j = ia[k]; j < ia[k+1] && !found; ++j)
				{
					if(ja[j] == rowi)
						a[j] += beta*elems[j0], found = true;
					else
						if(ja[j] < rowi) pos = j+1;
				}
				if(!found)
				{
					a.insert(a.begin()+pos, beta*elems[j0]);
					ja.insert(ja.begin()+pos, rowi);
					for(int j = k+1; j < n+1; ++j)
						ia[j]++;
				}
				++j0;
			}
		}
	}
	RemoveZerosRow(rowi);
}

void CSRMatrix::RemoveZeros(int row, bool flip_format)
{
	bool format = !flip_format ? csr : !csr;
	if(format)
	{
		std::vector<double>::iterator abeg = a.begin();
		std::vector<int>::iterator jabeg = ja.begin();
		int zeros = 0;
		for(int j = ia[row]; j < ia[row+1]; )
		{
			if(fabs(a[j]) < eps)
			{
				a.erase(a.begin()+j);
				ja.erase(ja.begin()+j);
				zeros++;
				ia[row+1]--;
			}
			else
				++j;
		}
		if(zeros)
			for(int k = row+2; k < n+1; ++k)
				ia[k] -= zeros;
	}
	else
	{
		for(int k = 0; k < n; ++k)
			for(int j = ia[k]; j < ia[k+1]; )
				if(ja[j] == row && fabs(a[j]) < eps)
				{
					a.erase(a.begin()+j);
					ja.erase(ja.begin()+j);
					for(int i = k+1; i < n+1; ++i)
						ia[i]--;
				}
				else
					++j;
	}

}
void CSRMatrix::RemoveZerosRow(int row)
{
	RemoveZeros(row, false);
}
void CSRMatrix::RemoveZerosCol(int col)
{
	RemoveZeros(col, true);
}
void CSRMatrix::RemoveZeros()
{
	for(int row = 0; row < n; ++row)
		RemoveZeros(row, !csr);
}
void CSRMatrix::AddRow(double alpha, int rowi, double beta, int rowj, bool in_place)
{
	Add(alpha, rowi, beta, rowj, in_place, false);
}
void CSRMatrix::AddCol(double alpha, int coli, double beta, int colj, bool in_place)
{
	Add(alpha, coli, beta, colj, in_place, true);
}
void CSRMatrix::PushElement(int row, int col, double elem)
{
	if(!csr)
	{
		int tmp = row; row = col; col = tmp;
	}
	int pos = ia[row];
	for(int j = ia[row]; j < ia[row+1]; ++j)
		if(ja[j] == col)
		{
			a[j] = elem;
			return;
		}
		else
			if(ja[pos] < col)
				pos = j+1;
	a.insert(a.begin()+pos, elem);
	ja.insert(ja.begin()+pos, col);
	for(int k = row+1; k < n+1; ++k)
		ia[k]++;
}


void CSRMatrix::ExtractDiagonal(std::vector<double>& diag) const
{
	diag.resize(n, 0.0);
	for(int i = 0; i < n; i++)
	{
		bool found = false;
		for(int k = ia[i]; !found && k < ia[i+1]; ++k)
			if(i == ja[k])
				diag[i] = a[k], found = true;
	}
}

double CSRMatrix::Diagonal(int i) const
{
	assert(i >= 0 && i < n);
	for(int j = ia[i]; j < ia[i+1]; ++j)
		if(ja[j] == i)
			return a[j];
	std::cerr << "Failed to find diagonal element" << std::endl;
	assert(false);
	return 0.0;
}

double CSRMatrix::L2Norm(int i) const
{
	double norm = 0.0;
	for(int j = ia[i]; j < ia[i+1]; ++j)
		norm += a[ja[j]]*a[ja[j]];
	return sqrt(norm);
}

double CSRMatrix::L1Norm(int i) const
{
	double norm = 0.0;
	for(int j = ia[i]; j < ia[i+1]; ++j)
		norm += fabs(a[ja[j]]);
	return norm;
}

bool CSRMatrix::Save(const char * filename, MatrixFormat fmt) const
{
	std::ofstream ofs(filename);
	int nnz = ia[n] - ia[0];
	ofs << n << " " << nnz << std::endl;
	for(int i = 0; i < n+1; ++i)
		ofs << (ia[i]+1) << " ";
	ofs << std::endl;
	for(int j = 0; j < nnz; ++j)
		ofs << (ja[j]+1) << " ";
	ofs << std::endl;
	for(int j = 0; j < nnz; ++j)
		ofs << a[j] << " ";
	ofs << std::endl;
	ofs.close();
	return true;
}


double CSRMatrix::Get(int row, int col) const
{
	for(int j = ia[row]; j < ia[row+1]; ++j)
		if(ja[j] == col) return a[j];
	return 0.0;
}

void RowAccumulator::SparseAdd(const std::vector<double>& a, const std::vector<int>& ja, double alpha, int jbeg, int jend)
{
        for(int j = jbeg; j < jend; ++j)
		Add(ja[j], alpha*a[j]);
}

void RowAccumulator::SparseAdd(const CSRMatrix* pA, int row, double alpha)
{
	SparseAdd(pA->GetA(), pA->GetJA(), alpha, pA->GetIA(row), pA->GetIA(row+1));
}

void RowAccumulator::SetRowFrom(const CSRMatrix* pA, int i, int beg, int end)
{
        int n = pA->Size();
        assert(i >= 0 && i < n);
	jr.resize(n,-1);
	int ia0 = pA->GetIA(i), ia1 = pA->GetIA(i+1);
	int nnz = ia0;
        for(int j = ia0; j < ia1; ++j)
        {
		int col = pA->GetJA(j);
		if(col >= beg && col < end)
		{
			jw.push_back(col);
			jr[col] = nnz-ia0;
			w.push_back(pA->GetA(j));
			++nnz;
		}
        }
}


void RowAccumulator::SetIntervalFrom(int n, const std::vector<int>& rows, const std::vector<double>& vals)
{
	if(jr.size() != n) jr.resize(n,-1);
	for(int i = 0; i < rows.size(); ++i)
	{
		int row = rows[i];
		jr[row] = jw.size();
		jw.push_back(row);
		w.push_back(vals[i]);
	}
}


void RowAccumulator::split(int p, int i)
{
	std::vector<int> perm(jw.size());
	for(int k = 0; k < perm.size(); ++k) perm[k] = k;
	split_rec(p, perm, 0);
}
void RowAccumulator::split_rec(int p, std::vector<int>& perm, int offset)
{
        int mid = partition(perm, offset);
	bool cycle = mid == offset && mid == p-1;
        if(mid >= p)
        {
                for(int k = mid; k < perm.size(); ++k)
			w[perm[k]] = 0.0;
		if(mid == p) return;
		else split_rec(p, perm, offset);
        }
	else split_rec(p, perm, mid+1);
}

// quick-sort algorithm
int RowAccumulator::partition(std::vector<int>& perm, int offset)
{
        int ipiv = offset;
        double piv = fabs(w[perm[offset]]);
        for(int k = perm.size()-1-offset; offset+k >= ipiv; )
        {
                if(fabs(w[perm[offset+k]]) > piv)
                {
                        int next = perm[offset+k];
                        perm[offset+k] = perm[ipiv+1];
                        for(int j = ipiv; j >= 0; --j)
                                perm[j+1] = perm[j];
                        perm[offset] = next;
                        ++ipiv;
                }
                else --k;
        }
        return ipiv;
}


void RowAccumulator::SelectLargest(int i, int p)
{
	if(w.size() <= p+1) return;
	split(p, i);
	RemoveZeros();
}

void RowAccumulator::Print() const
{
	for(int k = 0; k < jw.size(); ++k)
		std::cout << "col " << jw[k] << " val " << w[k]  << std::endl;
}

double RowAccumulator::L2Norm() const
{
	double norm = 0.0;
	for(int k = 0; k < w.size(); ++k)
		norm += w[k]*w[k];
	return sqrt(norm);
}

double RowAccumulator::L1Norm() const
{
	double norm = 0.0;
	for(int k = 0; k < w.size(); ++k)
		norm += fabs(w[k]);
	return norm;
}

int RowAccumulator::find_pos(int col) const
{
	for(int pos = 0; pos < jw.size(); ++pos)
		if(jw[pos] == col)
			return pos;
	return -1;
}

void RowAccumulator::Remove(int col)
{
	int pos = find_pos(col);
	if(pos != -1)
		remove(pos);
}

bool RowAccumulator::Drop(int col, double norm)
{
	int pos = find_pos(col);
	if(pos == -1)
		return false;
	else
	{
		assert(pos >= 0 && pos < w.size());
		if(fabs(w[pos]) < norm)
		{
			remove(pos);
			return true;
		}
		else return false;
	}
}

void RowAccumulator::Drop(int row, double norm, int p)
{
	for(int pos = jw.size()-1; pos >= 0; --pos)
		if(fabs(w[pos]) < norm && jw[pos] != row)
			remove(pos);
	SelectLargest(row, p);
}

void RowAccumulator::remove(int pos)
{
	for(int k = pos+1; k < jw.size(); ++k)
		--jr[jw[k]];
	jr[jw[pos]] = -1;
	jw.erase(jw.begin()+pos);
	w.erase(w.begin()+pos);
}

void RowAccumulator::RemoveZeros(double eps)
{
	for(int k = w.size()-1; k >= 0; --k)
		if(fabs(w[k]) < eps)
			remove(k);
}

double RowAccumulator::Get(int col) const
{
	if(jr[col] == -1)
		return 0.0;
	int pos = 0;
	while(pos < jw.size() && jw[pos] < col) ++pos;
	if(pos < w.size())
		return w[pos];
	else return 0.0;
}

bool RowAccumulator::set(int col, double val, bool add)
{
	if(jr[col] == -1) // new element
	{
		int pos = 0;
		while(pos < jw.size() && jw[pos] < col) ++pos;
		jw.insert(jw.begin()+pos, col);
		w.insert(w.begin()+pos, val);
		jr[col] = pos;
		for(int k = pos+1; k < jw.size(); ++k)
			++jr[jw[k]];
		return true;
	}
	else if(add) w[jr[col]] += val;
	else w[jr[col]] = val;
	return false;
}

void RowAccumulator::Scale(double alpha)
{
	if(fabs(alpha-1.0) > 1.0e-12)
		for(int i = 0; i < w.size(); ++i)
			w[i] *= alpha;
}

void RowAccumulator::Clear()
{
	for(int k = 0; k < jw.size(); ++k)
		jr[jw[k]] = -1;
	jw.clear();
	w.clear();
}

bool RowAccumulator::prepare_jr()
{
	bool prev = jralloc;
	if(!jralloc)
	{
		jr.resize(n);
		std::fill(jr.begin(), jr.end(), -1);
		for(int k = 0; k < jw.size(); ++k)
			jr[jw[k]] = k;
		jralloc = true;
	}
	return prev;
}

bool RowAccumulator::clear_jr()
{
	bool prev = jralloc;
	if(jralloc)
	{
		jr.clear();
		jralloc = false;
	}
	return prev;
}
