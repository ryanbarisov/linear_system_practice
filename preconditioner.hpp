#ifndef _PRECONDITIONER_H
#define _PRECONDITIONER_H

#include "matrix.hpp"
#include "priority_queue.hpp"
#include <map>
#include <utility>
#include <algorithm>


// #define AMG_COARSE_SOLVE // solve coarsest grid instead of smoothing



class SolverParameters
{
private:
	static const bool debug = true;
	std::map<std::string, double>	real_parameters;
	std::map<std::string, int>		integer_parameters;
public:
	void SetRealParameter(std::string name, double value)
	{
		if(debug)	std::cerr << "Set real parameter " << name << ": " << value  << std::endl;
		real_parameters[name] = value;
	}
	std::pair<bool,double> GetRealParameter(std::string name) const
	{
		static const std::pair<bool,double> parameter_not_found = std::make_pair(false, 0.0);
		std::map<std::string,double>::const_iterator it = real_parameters.find(name);
		if(it != real_parameters.end())	
			return std::make_pair(true, it->second);
		else
			return parameter_not_found;
	}
	void SetIntegerParameter(std::string name, int value)
	{
		if(debug)	std::cerr << "Set integer parameter " << name << ": " << value  << std::endl;
		integer_parameters[name] = value;
	}
	std::pair<bool,int> GetIntegerParameter(std::string name) const
	{
		static const std::pair<bool,int> parameter_not_found = std::make_pair(false, 0);
		std::map<std::string,int>::const_iterator it = integer_parameters.find(name);
		if(it != integer_parameters.end())	
			return std::make_pair(true, it->second);
		else
			return parameter_not_found;
	}
};




enum class PreconditionerType
{
	NONE,
	ILU0,
	ILUC,
	AMG
};

class Preconditioner
{
protected:
	PreconditionerType mytype;
	const MTX_matrix* pA;
public:
	Preconditioner(PreconditionerType _mytype, const MTX_matrix* _pA) : mytype(_mytype), pA(_pA) {}
	PreconditionerType GetType() const {return mytype;}

	virtual ~Preconditioner() {}
	
	virtual bool SetupPreconditioner() = 0;	
	virtual bool PreconditionedSolve(const std::vector<double>& rhs, std::vector<double>& x) = 0;
};



class ILUC_Preconditioner : public Preconditioner
{
private:
	S_matrix * L;
	S_matrix * U;

	double tau;
	int lfil;
public:
	ILUC_Preconditioner(const MTX_matrix* _pA, SolverParameters parameters) 
		: Preconditioner(PreconditionerType::ILUC, _pA), L(NULL), U(NULL) 
		{
			tau = parameters.GetRealParameter("drop_tolerance").second;
			lfil = parameters.GetIntegerParameter("level_of_fill").second;
		}

	~ILUC_Preconditioner() 
	{
		if(L != NULL)	delete L;
		if(U != NULL)	delete U;
	}

	bool SetupPreconditioner()
	{
		int n = pA->Size();
		L = new S_matrix(n, STORED_BY_COLS);
		U = new S_matrix(n, STORED_BY_ROWS);

		// prepare rows and columns of A into separate structures
		sparse_row * Arows = new sparse_row[n];
		sparse_row * Acols = new sparse_row[n];
		entry e;
		for(int k = 0; k < pA->L; k++)
		{
			int i = pA->ia[k]-1;
			int j = pA->ja[k]-1;
			e.second = pA->a[k];
			e.first = j;
			Arows[i].add_element(e);
			e.first = i;
			Acols[j].add_element(e);
		}
		// initialize bi-index structure
		int * Ufirst = new int[n];
		int * Lfirst = new int[n];
		std::list<entry> * Ulist = new std::list<entry>[n];
		std::list<entry> * Llist = new std::list<entry>[n];
		for(int k = 0; k < n; k++)
		{
			Ufirst[k] = Lfirst[k] = 0;
		}
		// while(Arows[0].row[Ufirst[0]].first < 1 && Ufirst[0] < Arows[0].row.size())	Ufirst[0]++;
		// while(Acols[0].row[Lfirst[0]].first < 1 && Lfirst[0] < Acols[0].row.size())	Lfirst[0]++;


		double * big_elems = new double[lfil];
		sparse_row z, w;
		for(int k = 0; k < n; k++)
		{
			z.assign(Arows[k],k,n);
			for(std::list<entry>::const_iterator it = Llist[k].begin(); it != Llist[k].end(); ++it)
				z.plus(-it->second, U->get_vector(it->first),k,n);
			w.assign(Acols[k],k+1,n);
			for(std::list<entry>::const_iterator it = Ulist[k].begin(); it != Ulist[k].end(); ++it)
				w.plus(-it->second, L->get_vector(it->first),k+1,n);
			
			double ukk = z.get_element(k);
			w *= 1.0/ukk;
			// TODO dropping rule for z and w

			// drop z
			for(int kk = 0; kk < lfil; kk++)
				big_elems[kk] = 0.0;
			for(sparse_type::iterator it = z.row.begin(); it != z.row.end(); )
			{
				double val = fabs(it->second);
				if(val < tau)
				{
					if(it->first == k)
					{
						std::cerr << "Delete diagonal element value " << val << " row " << k << std::endl;
						std::cin.ignore();
					}
					it = z.row.erase(it);
				}
				else
				{
					int kk = 0;
					while(kk < lfil && val < big_elems[kk])	kk++;
					if(kk < lfil)
					{
						for(int kkk = lfil-1; kkk > kk; kkk--)
							big_elems[kkk] = big_elems[kkk-1];
						big_elems[kk] = val;
					}
					it++;
				}
			}
			for(sparse_type::iterator it = z.row.begin(); it != z.row.end(); )
			{
				if(fabs(it->second) < big_elems[lfil-1] && it->first != k)
					it = z.row.erase(it);
				else it++;
			}
			// drop w
			for(int kk = 0; kk < lfil; kk++)
				big_elems[kk] = 0.0;
			for(sparse_type::iterator it = w.row.begin(); it != w.row.end(); )
			{
				double val = fabs(it->second);
				if(val < tau)
				{
					it = w.row.erase(it);
				}
				else
				{
					int kk = 0;
					while(kk < lfil && val < big_elems[kk])	kk++;
					if(kk < lfil)
					{
						for(int kkk = lfil-1; kkk > kk; kkk--)
							big_elems[kkk] = big_elems[kkk-1];
						big_elems[kk] = val;
					}
					it++;
				}
			}
			for(sparse_type::iterator it = w.row.begin(); it != w.row.end(); )
			{
				if(fabs(it->second) < big_elems[lfil-1])
					it = w.row.erase(it);
				else it++;
			}

			// update factors
			U->set_vector(k,z,k,n);
			L->set_vector(k,w,k+1,n);
			L->add_element(k,k,1.0);

			// update indices
			for(int i = 0; i < k+1; i++)
			{
				const sparse_row& rU = U->get_vector(i);
				if(rU.row[Ufirst[i]].first < k+1)
				{
					while(rU.row[Ufirst[i]].first < k+1 && Ufirst[i] < rU.row.size() ) Ufirst[i]++;
					if(Ufirst[i] < rU.row.size())
					{
						e.first = i;
						e.second = rU.row[Ufirst[i]].second;
						Ulist[ rU.row[Ufirst[i]].first ].push_back(e);
					}
				}
				const sparse_row& rL = L->get_vector(i);
				if(rL.row[Lfirst[i]].first < k+1)
				{
					while(rL.row[Lfirst[i]].first < k+1 && Lfirst[i] < rL.row.size() )	Lfirst[i]++;
					if(Lfirst[i] < rL.row.size())
					{
						e.first = i;
						e.second = rL.row[Lfirst[i]].second;
						Llist[ rL.row[Lfirst[i]].first ].push_back(e);
					}
				}
			}
		}

		delete [] big_elems;
		delete [] Ufirst;
		delete [] Lfirst;
		delete [] Llist;
		delete [] Ulist;
		delete [] Arows;
		delete [] Acols;

		// std::cout << "L:" << std::endl;
		// L->print();
		// std::cout << "U:" << std::endl;
		// U->print();

		// transpose L so it is faster to solve using LU
		S_matrix * Lt = new S_matrix(n, STORED_BY_ROWS), * L_tmp;
		for(int i = 0; i < n; i++)
		{
			const sparse_row& r = L->get_vector(i);
			for(sparse_type::const_iterator it = r.row.begin(); it != r.row.end(); ++it)
			{
				Lt->add_element(it->first, i, it->second);
			}
		}

		L_tmp = L; L = Lt; Lt = L_tmp;
		delete Lt;

		return true;
	}

	bool PreconditionedSolve(const std::vector<double>& rhs, std::vector<double>& x)
	{
		LU_solve(*L,*U,rhs,x);
		return true;
	}
};


////////////////////////////////////////////////////////////////////

class AMG_Preconditioner : public Preconditioner
{
	typedef std::map<int,sparse_row> interpolation_type;
private:
	int amg_levels = 2;

	std::vector<int> Nm;						// coarse system dimensions for 1 <= m <= q = amg_levels-1
	std::vector< std::vector<int> > Cm;			// index set for coarse grids, 1 <= m <= q = amg_levels-1
	std::vector<S_matrix*> Am;					// coarse grid matrices for 0 <= m <= amg_levels-1
	std::vector< interpolation_type > Wm;		// interpolation weights for 1 <= m <= amg_levels-1


	// prolongation from m to m-1, m > 0
	void prolongation(int m_from, const std::vector<double>& e, std::vector<double>& v_out)
	{
		assert(m_from > 0);
		int m = m_from - 1;
		int n = Am[m]->Size();
		v_out.resize(n);
		for(int i = 0; i < n; i++)	v_out[i] = 0.0;
		const std::vector<int>& C = Cm[m_from];
		interpolation_type& W = Wm[m_from];
		std::vector<bool> processed(n,false);
		std::map<int,int> invC;
		for(int k = 0; k < C.size(); k++)
		{
			int i = C[k];
			invC[i] = k;
			v_out[i] = e[k];
			processed[i] = true;
		}
		for(int i = 0; i < n; i++) if(!processed[i])
		{
			sparse_row& r = W[i];
			for(int k = 0; k < r.row.size(); k++)
			{
				int col = r.row[k].first;
				v_out[i] += r.row[k].second*e[invC[col]];
			}
		}
	}
	// restriction from m to m+1
	void restriction(int m_from, const std::vector<double>& e, std::vector<double>& v_out)
	{
		int m = m_from+1;
		int nf = e.size();
		int nc = Am[m]->Size();
		v_out.resize(nc);
		for(int i = 0; i < nc; i++)	v_out[i] = 0.0;
		const std::vector<int>& C = Cm[m];
		interpolation_type& W = Wm[m];
		std::vector<bool> processed(nf,false);
		// need to construct inverse mapping from C[k] to k to fill coarse system
		std::map<int,int> invC;
		for(int k = 0; k < C.size(); k++)
		{
			int j = C[k];
			invC[j] = k;
			v_out[k] = e[j];
			processed[j] = true;
		}
		for(int j = 0; j < nf; j++) if(!processed[j])
		{
			const sparse_row& r = W[j];
			for(int k = 0; k < r.row.size(); k++)
			{
				int row = r.row[k].first;
				v_out[invC[row]] += r.row[k].second*e[j];
			}
		}
	}

	void smoothing(int m, std::vector<double>& x, const std::vector<double>& b, int nrelax = 1)
	{
		int method = 0;	// 0 - Jacobi, 1 - Gauss-Seidel, 2 - ?
		for(int k = 0; k < nrelax; k++)
		{
			if(method == 0)
				jacobi(Am[m], x, b);
			else if(method == 1)
				gauss_seidel(Am[m], x, b);
		}
	}

	// construct system on step m+1 with given system on step m and interpolator matrices
	// A^{m+1} = I_m^{m+1} A^m I_{m+1}^m
	void construct_coarse_system(int m_from)
	{
		const S_matrix* Amf = Am[m_from];	// fine system

		int m = m_from+1;
		int nc = Nm[m];					// coarse system size
		int nf = Am[m_from]->Size();		// fine system size

		// (1) A^{m+1} = I_m^{m+1} A^m
		// (2) A^{m+1} = A^{m+1} I_{m+1}^m
		
		std::cout << "Fine system level " << m_from << " size " << Amf->Size() << std::endl;
		std::cout << "Coarse system level " << m << " size " << nc << std::endl;


		S_matrix* Amc = Am[m];
		Amc = new S_matrix(nc, true);	// matrix stored by rows

		// need to construct inverse mapping from C[k] to k to fill coarse system
		std::map<int,int> invC;

		// step (1)
		const std::vector<int>& C = Cm[m];
		assert(nc == C.size());
		interpolation_type& W = Wm[m];
		std::vector<bool> coarse(nf,false);
		// coarse indices 
		for(int k = 0; k < C.size(); k++)
		{
			int j = C[k];
			invC[j] = k;
			coarse[j] = true;
		}
		for(int k = 0; k < C.size(); k++)
		{
			int j = C[k];
			Amc->set_vector(k, Amf->get_vector(j), 0, nf);
			// sparse_row& r = Amc->get_vector(k);
			// std::cout << "coarse " << k << " vector: " << std::endl;
			// r.print();
			// std::cin.ignore();
		}
		// fine indices
		for(int j = 0; j < nf; j++) if(!coarse[j])
		{
			const sparse_row& r = W[j];
			// std::cout << "Weights row " << j << std::endl;
			// r.print();
			// assert(!r.row.empty());
			for(int k = 0; k < r.row.size(); k++)
			{
				int row = r.row[k].first;
				// std::cout << "coarse " << invC[row] << ":" << std::endl;
				// Amc->get_vector(invC[row]).print();
				// std::cout << " add vector " << row << " with coef " << r.row[k].second << ":" << std::endl;
				// Amf->get_vector(row).print();
				Amc->add_vector(invC[row], r.row[k].second, Amf->get_vector(j), 0, nf);

				// std::cout << "after adding coarse " << invC[row] << ":" << std::endl;
				// Amc->get_vector(invC[row]).print();
				// std::cin.ignore();
			}
		}

		// Amc->print();
		// std::cin.ignore();

		// step (2)
		S_matrix* A = new S_matrix(nc, true);	// temporary intermediate matrix 
		for(int Arow = 0; Arow < nc; Arow++)
		{
			const sparse_row& r = Amc->get_vector(Arow);
			sparse_row& r2 = A->get_vector(Arow);
			for(int k = 0; k < r.row.size(); k++)
			{
				int Acol = r.row[k].first;
				double Aelem = r.row[k].second;
				if(coarse[Acol])
				{
					A->add_element(Arow, invC[Acol], Aelem);
					// A->add_element(Arow, Acol, Aelem);
				}
				else
				{
					const sparse_row& rW = W[Acol];
					for(int kk = 0; kk < rW.row.size(); kk++)
					{
						int col = rW.row[kk].first;
						if(fabs(rW.row[kk].second) > 1e-10)
						{
							// std::cout << "A " << Arow << " col " << col 
							// 	<< " invC(col) " << invC[col] << " elem " << rW.row[kk].second * Aelem
							// 	<< std::endl;
							A->add_element(Arow, invC[col], rW.row[kk].second * Aelem);
						}
					}
				}
			}
			r2.remove_zeros();
		}

		S_matrix* tmpA = Amc;
		Amc = A;
		A = tmpA;

		delete tmpA;

		Am[m] = Amc;
		// Am[m]->print();

		char str[256];
		sprintf(str, "amg_level_%02d.mtx", m);
		Am[m]->Save(str);
	}

	void construct_coarsest_inverse()
	{
		S_matrix* A = Am[amg_levels-1];
		int n = A->Size();

		std::vector<double> diag(n);
		for(int k = 1; k < n; k++)
		{
			const sparse_row& r = A->get_vector(k);
			for(int kk = 0; kk < r.row.size(); kk++)
			{
				int i = r.row[kk].first;
				if(k == i)
				{
					diag[k] = r.row[kk].second;
					break;
				}
			}
		}
		// ilu0
		for(int i = 1; i < n; i++)	// loop over rows
		{
			double akk, aik;
			int kb = 0, ke = i;
			sparse_row& r = A->get_vector(i);
			for(int kk = 0; kk < r.row.size(); kk++)
			{
				int k = r.row[kk].first;
				if(!(k >= kb && k > ke))	continue;
				r.row[kk].second /= diag[k];
				aik = r.row[kk].second;

				const sparse_row& r2 = A->get_vector(k);
				int jb = k+1, je = n;
				for(int kkk = kk+1; kkk < r.row.size(); kkk++)
				{
					int j = r.row[kkk].first;
					if(!(j >= jb && j < je))	continue;
					for(int kkkk = 0; kkkk < r2.row.size(); kkkk++)
					{
						if(r2.row[kkkk].first == j)
						{
							r.row[kkk].second -= aik*r2.row[kkkk].second;
						}
					}
				}
			}
		}
	}

	void solve_coarsest(std::vector<double>& x, const std::vector<double>& b)
	{
		S_matrix* A = Am[amg_levels-1];
		int n = A->Size();

		// solve Lx = b
		for(int i = 0; i < n; i++)	// loop over rows
		{
			const sparse_row& r = A->get_vector(i);
			x[i] = b[i];
			for(int kk = 0; kk < r.row.size(); kk++)
			{
				int k = r.row[kk].first;
				if(k < i)
					x[i] -= r.row[kk].second*x[k];
			}
		}
		// solve Uy = x
		for(int i = n-1; i >= 0; i--)	// loop over rows
		{
			double akk = 0.0;
			const sparse_row& r = A->get_vector(i);
			for(int kk = 0; kk < r.row.size(); kk++)
			{
				int k = r.row[kk].first;
				if(k == i)
				{
					akk = r.row[kk].second;
					break;
				}
			}
			assert(akk != 0.0);
			x[i] = b[i];
			for(int kk = 0; kk < r.row.size(); kk++)
			{
				int k = r.row[kk].first;
				if(k > i)
					x[i] -= r.row[kk].second*x[k];
			}
			x[i] /= akk;
		}
	}

	void find_strongly_connected(const S_matrix* pA, std::vector<std::vector<int> >& sc, double theta = 0.25)
	{
		int n = pA->Size(), row, col, value;
		std::vector<double> maxmod(n, 0.0);
		sc.resize(n);

		for(int row = 0; row < n; row++)
		{
			const sparse_row& r = pA->v[row];
			for(int j = 0; j < r.row.size(); j++)
			{
				col = r.row[j].first;
				value = fabs(r.row[j].second);
				if(row != col && maxmod[row] < value)	maxmod[row] = value;
			}
		}

		for(int row = 0; row < n; row++)
		{
			const sparse_row& r = pA->v[row];
			std::vector<int>& sci = sc[row];
			double rowmax = maxmod[row];
			for(int j = 0; j < r.row.size(); j++)
			{
				col = r.row[j].first;
				value = fabs(r.row[j].second);
				if(value > theta * rowmax)	sci.push_back(col);
			}
			// sort to use binary search over strong connections later
			if(!sc[row].empty())	std::sort(sc[row].begin(), sc[row].end());
		}
	}

	void find_adjacent_strongly_connected(const std::vector<std::vector<int> >& sc, std::vector<std::vector<int> >& asc)
	{
		int n = sc.size(), m, i;
		for(int j = 0; j < n; j++)
		{
			const std::vector<int> v = sc[j];
			m = v.size();
			for(int k = 0; k < m; k++)
			{
				i = v[k];
				asc[i].push_back(j);
			}
		}
	}

public:
	AMG_Preconditioner(const MTX_matrix* _pA, SolverParameters parameters) 
		: Preconditioner(PreconditionerType::AMG, _pA)
		{
			amg_levels = std::max(2, parameters.GetIntegerParameter("amg_levels").second);
			Nm.resize(amg_levels);
			Cm.resize(amg_levels);
			Am.resize(amg_levels);
			Wm.resize(amg_levels);
			Am[0] = new S_matrix(*_pA);
			Nm[0] = Am[0]->Size();
			for(int i = 0; i < Nm[0]; i++)	Cm[0].push_back(i);	// is it needed?
		}

	AMG_Preconditioner(const AMG_Preconditioner& other) : Preconditioner(PreconditionerType::AMG, other.pA), 
		amg_levels(other.amg_levels), Nm(other.Nm), Cm(other.Cm), Wm(other.Wm)
	{
		Am.resize(amg_levels);
		//Am[0] = new S_matrix(*(other.pA));
		for(int i = 0; i < amg_levels; i++)
			Am[i] = new S_matrix( *(other.Am[i]) );
	}

	AMG_Preconditioner& operator=(const AMG_Preconditioner& other)
	{
		mytype = PreconditionerType::AMG;
		pA = other.pA;
		amg_levels = other.amg_levels;
		Nm = other.Nm;
		Cm = other.Cm;
		Wm = other.Wm;
		for(int i = 0; i < amg_levels; i++)
			if(Am[i] != NULL)	delete Am[i];
		Am.resize(amg_levels);
		// Am[0] = new S_matrix(*(other.pA));
		for(int i = 0; i < amg_levels; i++)
			Am[i] = new S_matrix( *(other.Am[i]) );
	}

	~AMG_Preconditioner() 
	{
		for(int i = 0; i < amg_levels; i++)
			if(Am[i] != NULL)	delete Am[i];
	}


	void construct_CF_partition(int m_from)
	{
		bool debug = true;

		int m_to = m_from + 1;
		int m = m_to;

		const S_matrix* pmat = Am[m_from];
		int n = pmat->Size();
		PriorityQueue<int, std::greater<int> > q(n);
		std::vector< std::vector<int> > strong_connections(n);		// strong connections for each grid element
		std::vector< std::vector<int> > adj_strong_connections(n);	// adjacent strong connections for each element

		std::vector<int>& C = Cm[m];	// coarse	partition on level m
		std::vector<int> F;				// fine 	partition on level m

		find_strongly_connected(pmat, strong_connections);
		find_adjacent_strongly_connected(strong_connections, adj_strong_connections);
		for(int i = 0; i < adj_strong_connections.size(); i++)
			q.Push(i, adj_strong_connections[i].size());

		int max_asc_amount = n;	// maximum amount of adjacent strong connections cannot be more than matrix size
		while(!q.Empty())
		{
			int i = q.Peek();
			C.push_back(i);
			q.Pop();
			const std::vector<int>& asc_i = adj_strong_connections[i];
			for(int k = 0; k < asc_i.size(); k++)
			{
				int j = asc_i[k];
				if(q.Contains(j))
				{
					F.push_back(j);
					// remove j from q
					q.ChangeKey(j, n);
					assert( q.Pop() == j );
					const std::vector<int>& sc_j = strong_connections[j];
					for(int kk = 0; kk < sc_j.size(); kk++)
					{
						int l = sc_j[kk];
						if(q.Contains(l))	
							q.ChangeKey(l, q.GetKey(l)+1);
							// q.IncreaseKey(l, q.GetKey(l)+1);
					}
				}
			}
			const std::vector<int>& sc_i = strong_connections[i];
			for(int k = 0; k < sc_i.size(); k++)
			{
				int j = sc_i[k];
				if(q.Contains(j))	
					q.ChangeKey(j, q.GetKey(j)-1);
					// q.DecreaseKey(j, q.GetKey(j)-1);
			}
		}
	
		// extend constructed coarse partition	
		std::sort(C.begin(), C.end());
		std::sort(F.begin(), F.end());
		std::vector<bool> F_available(F.size(), true);
		int F_amount = F.size();	// amount of fine mesh elements

		if(false)
		{
			std::cout << "before extension coarse " << C.size() << std::endl;
			for(int i = 0; i < C.size(); i++)
			{
				std::cout << C[i] << " ";
			}
			std::cout << std::endl;
			std::cout << "before extension fine " << F.size() << std::endl;
			for(int i = 0; i < F.size(); i++)
			{
				std::cout << F[i] << " ";
			}

			std::cout << std::endl;
			std::cin.ignore();
		}

		for(int k = 0; k < F.size(); k++)
		{
			if(!F_available[k])	continue;	// skip disabled elements 
			int i = F[k];
			std::vector<int> hat_Cm;
			std::vector<int> Im;	// interpolation connections
			std::vector<int> Ds;	// skipped strong connections
			const std::vector<int>& sc_i = strong_connections[i];
			// initialize connections
			for(int kk = 0; kk < sc_i.size(); kk++)
			{
				int l = sc_i[kk];
				if(std::binary_search(C.begin(), C.end(), l))
					Im.push_back(l);
				else
					Ds.push_back(l);
			}
			std::sort(Im.begin(), Im.end());

			for(int kk = 0; kk < Ds.size(); kk++)
			{
				int j = Ds[kk];
				const std::vector<int>& sc_j = strong_connections[j];
				// intersect sc_j with Im \cup hat_Cm
				std::vector<int> v_intersection;
				std::set_intersection(sc_j.begin(), sc_j.end(), Im.begin(), Im.end(), std::back_inserter(v_intersection));
				if(v_intersection.empty())
				{
					std::vector<int>::iterator jt = Im.begin();
					while(jt != Im.end() && j >= *jt)	++jt;
					Im.insert(jt, j);
					hat_Cm.push_back(j);
				}
			}
			if(hat_Cm.size() > 1)
			{
				std::vector<int>::iterator jt = C.begin();
				while(jt != C.end() && i >= *jt)	++jt;
				C.insert(jt, i);
				F_available[k] = false;
				F_amount--;
			}
			else if(!hat_Cm.empty())
			{
				int l = hat_Cm[0];
				std::vector<int>::iterator jt = C.begin();
				while(jt != C.end() && l >= *jt)	++jt;
				C.insert(jt, l);

				int kk = 0;
				while(kk < F.size() && l < F[kk])	++kk;
				if(kk < F.size() && l == F[kk])
				{
					F_available[kk] = false;
					F_amount--;
				}
			}
		}

		if(false)
		{
			std::cout << "after extension coarse " << C.size() << " fine " << F_amount << std::endl;
			for(int i = 0; i < C.size(); i++)
			{
				std::cout << C[i] << " ";
			}
			std::cout << std::endl;
			std::cin.ignore();
		}

		Nm[m] = C.size();

		// repartition F so that disabled elements appear at the end
		int nF = F.size(), itmp;
		int index1 = 0, index2 = nF-1;
		while(index1 < index2)
		{
			while(index1 < index2 &&  F_available[index1])	index1++;
			while(index1 < index2 && !F_available[index2])	index2--;
			if(index1 < index2)
			{
				itmp = F[index1];
				F[index1] = F[index2];
				F[index2] = itmp;
				F_available[index1] = true;
				F_available[index2] = false;
			}
			else if(index1 == nF-1)	index1++;
		}
		// remove disabled elements and sort 
		nF = index1;
		F.resize(nF);
		std::sort(F.begin(), F.end());

		if(false)
		{
			std::cout << "fine " << F.size() << std::endl;
		}

		// determine interpolation weights
		interpolation_type& W = Wm[m];

		int col;
		double val, di, sj, dk;
		for(int k = 0; k < F.size(); k++)
		{
			int i = F[k];
			// std::map<int,double> Im;	// interpolation connections and intermediate coefficients
			std::vector<int> Im_col;
			std::vector<double> Im_coef;
			std::vector<int> Ds;	// skipped strong connections
			// std::vector<int> Dw;	// skipped weak connections
			const std::vector<int>& sc_i = strong_connections[i];
			// initialize connections
			// strong_connections is sorted, that's why Im and Ds are also sorted 
			const sparse_row& ri = pmat->v[i];
			for(int kk = 0; kk < sc_i.size(); kk++)
			{
				int l = sc_i[kk];
				if(std::binary_search(C.begin(), C.end(), l))
				{
					// find element with column j
					val = 0.0;
					for(int kkk = 0; kkk < ri.row.size(); kkk++)
						if(ri.row[kkk].first == l)
						{
							val = ri.row[kkk].second;
							break;
						}
					if(val != 0.0)
					{
						Im_col.push_back(l);
						Im_coef.push_back(val);
					}
				}
				else
					Ds.push_back(l);
			}
			for(int kk = 0; kk < Im_col.size(); kk++)
				for(int kkk = kk+1; kkk < Im_col.size(); kkk++)
				{
					if(Im_col[kk] > Im_col[kkk])
					{
						col = Im_col[kk];
						Im_col[kk] = Im_col[kkk];
						Im_col[kkk] = col;
						val = Im_coef[kk];
						Im_coef[kk] = Im_coef[kkk];
						Im_coef[kkk] = val;
					}
				}
			
			if(false)
			{
				std::cout << "Im " << i << ":\n";
				for(int kk = 0; kk < Im_col.size(); kk++)
				{
					std::cout << " col " << Im_col[kk] << " val " << Im_coef[kk];
				}
				std::cout << std::endl;
				std::cout << "Ds " << i << ":\n";
				for(int kk = 0; kk < Ds.size(); kk++)
				{
					std::cout << " " << Ds[kk];
				}
				std::cout << std::endl;
				std::cin.ignore();
			}
			
			di = 0.0;
			for(int kk = 0; kk < ri.row.size(); kk++)
			{
				col = ri.row[kk].first;
				val = ri.row[kk].second;
				if(col == i)
					di += val;
				else if(!std::binary_search(sc_i.begin(), sc_i.end(), col))
				{
					// Dw.push_back(col);
					di += val;
				}
			}
			for(int kk = 0; kk < Ds.size(); kk++)
			{
				int j = Ds[kk];
				const std::vector<int>& sc_j = strong_connections[j];
				const sparse_row& rj = pmat->v[j];
				sj = 0.0;
				// std::vector<std::pair<int,double> > columns_values_Im_cap_Sj;
				std::map<int,double> columns_values_Im_cap_Sj;
				std::vector<int> intersection;	// work vector to store sorted array intersections
				std::set_intersection(sc_j.begin(),sc_j.end(),Im_col.begin(),Im_col.end(),std::back_inserter(intersection));
				for(int kkk = 0; kkk < rj.row.size(); kkk++) 
				{
					col = rj.row[kkk].first;
					if(std::binary_search(intersection.begin(),intersection.end(),col))	
					{
						double ajk = rj.row[kkk].second;
						sj += ajk;
						val = 0.0;
						for(int kkkk = 0; kkkk < ri.row.size(); kkkk++)
							if(ri.row[kkkk].first == j)
							{
								val = ri.row[kkkk].second;
								break;
							}
						// columns_values_Im_cap_Sj.push_back(std::make_pair(col, val*ajk));
						// std::cout << "update column " << col << " j " << j << " k " << col << " aij " << val << " ajk " << ajk << std::endl;
						columns_values_Im_cap_Sj[col] += val*ajk;
					}
				}
				if(sj != 0.0)
				{
					std::vector<int>::const_iterator it_beg = Im_col.begin(), it_end = Im_col.end();
					// for(int kkk = 0; kkk < columns_values_Im_cap_Sj.size(); kkk++)
					for(std::map<int,double>::const_iterator itt = columns_values_Im_cap_Sj.begin(); itt != columns_values_Im_cap_Sj.end(); itt++)
					{
						std::vector<int>::const_iterator it = std::find(it_beg,it_end,itt->first);
						assert(it != it_end);
						int ind = it-it_beg;

						// std::cout << "update " << ind << " before " << Im_coef[ind] << " += " << itt->second << "/" << sj << std::endl;
						Im_coef[ind] += itt->second / sj;
						// std::cout << "update " << ind << " after " << Im_coef[ind] << std::endl;
					}
				}
				// std::cout << "Ds " << i << " index " << j << " sj " << sj << std::endl;
			}
			if(!Im_col.empty())
			{
				sparse_row& rW = W[i];
				for(int kk = 0; kk < Im_col.size(); kk++)
				{
					rW.add_element(std::make_pair(Im_col[kk], -Im_coef[kk] / di));
					// rW.add_element(std::make_pair(Im_col[kk], 1.0));
					// std::cout << "i " << i << " col " << Im_col[kk] << " val " 
					// 	<< -Im_coef[kk] / di << " di " << di << " it->second " << Im_coef[kk] << std::endl;
				}
			}
			else
			{
				//?
			}
			// std::cin.ignore();
		}
	}
	



	bool SetupAMG(int m)
	{
		bool stop = m >= amg_levels-1;
		if(!stop)
		{
			construct_CF_partition(m);
			construct_coarse_system(m);
			SetupAMG(m+1);
		}
#if defined(AMG_COARSE_SOLVE)
		else
			construct_coarsest_inverse();
#endif
		return true;
	}

	bool SetupPreconditioner()
	{
		return SetupAMG(0);
	}

private:
	void V_cycle(int m, std::vector<double>& x, const std::vector<double>& b)
	{
		const S_matrix* A = Am[m];
		bool stop = m >= amg_levels-2;
		int niters_smooth = 2;
		// pre-smoothing
		smoothing(m,x,b, niters_smooth);
		// compute residual r = b - Ax
		int nc = Am[m+1]->Size();
		int n = A->Size();
		std::vector<double> resid(n), resid_restricted(nc), x0(nc, 0.0);
		for(int i = 0; i < n; i++)
		{
			const sparse_row& r = A->get_vector(i);
			resid[i] = b[i];
			for(int j = 0; j < r.row.size(); j++)
				resid[i] -= r.row[j].second * x[r.row[j].first];
		}
		restriction(m,resid,resid_restricted);
		if(!stop)
			V_cycle(m+1,x0,resid_restricted);
		else
		{
#if defined(AMG_COARSE_SOLVE)
			solve_coarsest(x0,resid_restricted);
#else
			smoothing(m+1,x0,resid_restricted, 1);
#endif
		}	
		prolongation(m+1,x0,resid);
		for(int i = 0; i < n; i++)
			x[i] += resid[i];
		// post-smoothing
		smoothing(m,x,b, niters_smooth);
	}

	void W_cycle(int m, std::vector<double>& x, const std::vector<double>& b)
	{
		const S_matrix* A = Am[m];
		bool stop = m >= amg_levels-2;
		int niters_smooth = 2;
		// pre-smoothing
		smoothing(m,x,b, niters_smooth);
		// compute residual r = b - Ax
		int nc = Am[m+1]->Size();
		int n = A->Size();
		std::vector<double> resid(n), resid_restricted(nc), x0(nc, 0.0);
		for(int i = 0; i < n; i++)
		{
			const sparse_row& r = A->get_vector(i);
			resid[i] = b[i];
			for(int j = 0; j < r.row.size(); j++)
				resid[i] -= r.row[j].second * x[r.row[j].first];
		}
		restriction(m,resid,resid_restricted);
		if(!stop)
			W_cycle(m+1,x0,resid_restricted);
		else
		{
#if defined(AMG_COARSE_SOLVE)
			solve_coarsest(x0,resid_restricted);
#else
			smoothing(m+1,x0,resid_restricted, 1);
#endif
		}
		prolongation(m+1,x0,resid);
		for(int i = 0; i < n; i++)
			x[i] += resid[i];
		// re-smoothing
		smoothing(m,x,b, niters_smooth);
		// compute residual r = b - Ax
		for(int i = 0; i < n; i++)
		{
			const sparse_row& r = A->get_vector(i);
			resid[i] = b[i];
			for(int j = 0; j < r.row.size(); j++)
				resid[i] -= r.row[j].second * x[r.row[j].first];
		}
		restriction(m,resid,resid_restricted);
		if(!stop)
			W_cycle(m+1,x0,resid_restricted);
		else
		{
#if defined(AMG_COARSE_SOLVE)
			solve_coarsest(x0,resid_restricted);
#else
			smoothing(m+1,x0,resid_restricted, 1);
#endif
		}
		prolongation(m+1,x0,resid);
		for(int i = 0; i < n; i++)
			x[i] += resid[i];
		// post-smoothing
		smoothing(m,x,b, niters_smooth);
	}
public:

	bool PreconditionedSolve(const std::vector<double>& rhs, std::vector<double>& x)
	{
		int n = x.size();
		for(int i = 0; i < n; i++)	x[i] = 0.0;
		W_cycle(0, x, rhs);	
		return true;
	}
};





////////////////////////////////////////////////////////////////////









static Preconditioner * CreatePreconditioner(PreconditionerType ptype, const MTX_matrix* pA, const SolverParameters& parameters)
{
	if(ptype == PreconditionerType::ILUC)
	{
		return new ILUC_Preconditioner(pA, parameters);
	}
	else if(ptype == PreconditionerType::AMG)
	{
		return new AMG_Preconditioner(pA, parameters);
	}
	else
	{
		std::cerr << "Sorry, preconditioner was not implemented " << std::endl;
		return NULL;
	}
}


#endif