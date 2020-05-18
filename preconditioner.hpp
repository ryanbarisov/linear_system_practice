#ifndef _PRECONDITIONER_H
#define _PRECONDITIONER_H

#include "matrix.hpp"
#include "priority_queue.hpp"
#include <map>
#include <utility>
#include <algorithm>

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

// forward declarations
class ILUC_Preconditioner;	


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
	void prolongation(int _m, const std::vector<double>& e, std::vector<double>& v_out)
	{
		assert(_m > 0);
		int m = _m - 1;
		int n = Am[m]->Size();
		v_out.resize(n);
		const std::vector<int>& C = Cm[m];
		interpolation_type& W = Wm[m];
		std::vector<bool> processed(n,false);
		for(int k = 0; k < C.size(); k++)
		{
			int i = C[k];
			v_out[i] = e[k];
			processed[i] = true;
		}
		for(int i = 0; i < n; i++) if(!processed[i])
		{
			sparse_row& r = W[i];
			assert(!r.row.empty());
			v_out[i] = r.sparse_dot(e);
		}
	}
	// restriction from m to m+1
	void restriction(int _m, const std::vector<double>& e, std::vector<double>& v_out)
	{
		int m = _m;
		int n = Am[m+1]->Size();
		v_out.resize(n,0.0);
		const std::vector<int>& C = Cm[m];
		interpolation_type& W = Wm[m];
		std::vector<bool> processed(n,false);
		for(int k = 0; k < C.size(); k++)
		{
			int j = C[k];
			v_out[j] = e[k];
			processed[j] = true;
		}
		for(int j = 0; j < n; j++) if(!processed[j])
		{
			const sparse_row& r = W[j];
			assert(!r.row.empty());
			for(int k = 0; k < r.row.size(); k++)
			{
				int row = r.row[k].first;
				v_out[row] += r.row[k].second*e[j];
			}
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

	void initialize_priority_queue(const std::vector<std::vector<int> >& asc, PriorityQueue<int>& q)
	{
		int n = asc.size();
		for(int i = 0; i < n; i++)
			q.Push(i, asc[i].size());
	}

public:
	AMG_Preconditioner(const MTX_matrix* _pA, SolverParameters parameters) 
		: Preconditioner(PreconditionerType::AMG, _pA)
		{
			amg_levels = std::max(2, parameters.GetIntegerParameter("amg_levels").second);
			Nm.resize(amg_levels-1);
			Cm.resize(amg_levels-1);
			Am.resize(amg_levels);
			Wm.resize(amg_levels-1);
			Am[0] = new S_matrix(*_pA);
		}

	AMG_Preconditioner(const AMG_Preconditioner& other) : Preconditioner(PreconditionerType::AMG, other.pA), 
		amg_levels(other.amg_levels), Nm(other.Nm), Cm(other.Cm), Wm(other.Wm)
	{
		Am.resize(amg_levels);
		Am[0] = new S_matrix(*(other.pA));
		for(int i = 1; i < amg_levels; i++)
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
		Am[0] = new S_matrix(*(other.pA));
		for(int i = 1; i < amg_levels; i++)
			Am[i] = new S_matrix( *(other.Am[i]) );
	}

	~AMG_Preconditioner() 
	{
		for(int i = 0; i < amg_levels; i++)
			if(Am[i] != NULL)	delete Am[i];
	}


	void construct_CF_partition(int m)
	{
		const S_matrix* pmat = Am[m];
		int n = pmat->Size();
		PriorityQueue<int> q(n);
		std::vector< std::vector<int> > strong_connections(n);		// strong connections for each grid element
		std::vector< std::vector<int> > adj_strong_connections(n);	// adjacent strong connections for each element

		std::vector<int>& C = Cm[m];	// coarse	partition on level m
		std::vector<int> F;				// fine 	partition on level m

		find_strongly_connected(pmat, strong_connections);
		find_adjacent_strongly_connected(strong_connections, adj_strong_connections);
		initialize_priority_queue(adj_strong_connections, q);

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
					q.IncreaseKey(j, max_asc_amount);
					q.Pop();
					const std::vector<int>& sc_j = strong_connections[j];
					for(int kk = 0; kk < sc_j.size(); kk++)
					{
						int l = sc_j[kk];
						if(q.Contains(l))	q.IncreaseKey(l, q.GetKey(l)+1);
					}
				}
			}
			const std::vector<int>& sc_i = strong_connections[i];
			for(int k = 0; k < sc_i.size(); k++)
			{
				int j = sc_i[k];
				if(q.Contains(j))	q.DecreaseKey(j, q.GetKey(j)-1);
			}
		}
	
		// extend constructed coarse partition	
		std::sort(C.begin(), C.end());
		std::sort(F.begin(), F.end());
		std::vector<bool> F_available(F.size(), true);
		int F_amount = F.size();	// amount of find mesh elements

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
			else
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

		Nm[m] = C.size();

		// repartition F so that disabled elements to appear at the end
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

		// determine interpolation weights
		interpolation_type& W = Wm[m];

		int col;
		double val, di, sj, dk;
		for(int k = 0; k < F.size(); k++)
		{
			int i = F[k];
			std::map<int,double> Im;	// interpolation connections and intermediate coefficients
			std::vector<int> Ds;	// skipped strong connections
			std::vector<int> Dw;	// skipped weak connections
			const std::vector<int>& sc_i = strong_connections[i];
			// initialize connections
			// strong_connections is sorted, that's why Im and Ds are also sorted 
			for(int kk = 0; kk < sc_i.size(); kk++)
			{
				int l = sc_i[kk];
				if(std::binary_search(C.begin(), C.end(), l))
				{
					// find element with column j
					val = 0.0;
					const sparse_row& r = pmat->v[i];
					for(int kkk = 0; kkk < r.row.size() && val != 0.0; kkk++)
						if(r.row[kkk].first == l)	val = r.row[kkk].second;
					Im.insert(std::make_pair(l,val));
				}
				else
					Ds.push_back(l);
			}
			
			di = 0.0;
			const sparse_row& ri = pmat->v[i];
			for(int kk = 0; kk < ri.row.size(); kk++)
			{
				col = ri.row[kk].first;
				val = ri.row[kk].second;
				if(col == i)
					di += val;
				else if(fabs(val) > 1e-15 && col != i && !std::binary_search(sc_i.begin(), sc_i.end(), col))
				{
					Dw.push_back(col);
					di += val;
				}
			}
			for(int kk = 0; kk < Ds.size(); kk++)
			{
				int j = Ds[kk];
				const std::vector<int>& sc_j = strong_connections[j];
				const sparse_row& rj = pmat->v[j];
				sj = 0.0;
				std::vector<std::pair<int,double> > columns_values_Im_cap_Sj;
				for(int kkk = 0; kkk < rj.row.size(); kkk++) 
				{
					col = rj.row[kkk].first;
					if(Im.find(col) != Im.end() && std::binary_search(sc_j.begin(),sc_j.end(),col))	
					{
						val = rj.row[kkk].second;
						sj += val;
						val = 0.0;
						for(int kkkk = 0; kkkk < ri.row.size(); kkkk++)
							if(ri.row[kkkk].first == j)	val = ri.row[kkkk].second;
						columns_values_Im_cap_Sj.push_back(std::make_pair(col, val*rj.row[kkk].second));
					}
				}
				if(sj != 0.0)
				{
					for(int kkk = 0; kkk < columns_values_Im_cap_Sj.size(); kkk++)
						Im[columns_values_Im_cap_Sj[kkk].first] += columns_values_Im_cap_Sj[kkk].second / sj;
				}
			}
			sparse_row& rW = W[i];
			for(std::map<int,double>::const_iterator it = Im.begin(); it != Im.end(); ++it)
			{
				rW.add_element(std::make_pair(it->first, -it->second / di));
			}
		}
	}
	



	bool SetupAMG(int m)
	{
		construct_CF_partition(m);
		// construct_interpolator(m);
		// construct_system(m+1);
		bool stop = m < amg_levels-1;
		if(!stop)
		{
			SetupAMG(m+1);
		}
		else
		{
			// construct_inverse();
		}
		return true;
	}

	bool SetupPreconditioner()
	{
		return SetupAMG(0);
	}

	bool PreconditionedSolve(const std::vector<double>& rhs, std::vector<double>& x)
	{
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
	else 
	{
		std::cerr << "Sorry, preconditioner was not implemented " << std::endl;
		return NULL;
	}
}


#endif