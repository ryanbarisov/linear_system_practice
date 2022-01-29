#include <preconditioner.h>
#include <solver_parameters.hpp>
#include <matrix.h>
#include <method.h>
#include <algorithm>
#include <set>

//#define AMG_COARSE_SOLVE // solve coarsest grid instead of smoothing


AMG_Preconditioner::AMG_Preconditioner(const SparseMatrix* pA, SolverParameters parameters) 
	: Preconditioner(PreconditionerType::AMG, pA)
	{
		amg_levels = std::max(2, parameters.GetIntegerParameter("amg_levels").second);
		niters_smooth = std::max(0, parameters.GetIntegerParameter("niters_smooth").second);
		w_cycle = std::max(0, parameters.GetIntegerParameter("w_cycle").second);
		Nm.resize(amg_levels);
		Cm.resize(amg_levels);
		Am.resize(amg_levels);
		Wm.resize(amg_levels);
		Am[0] = new SparseMatrix(*pA);
		Nm[0] = Am[0]->Size();
		for(int i = 0; i < Nm[0]; i++)	Cm[0].push_back(i);	// is it needed?
	}

AMG_Preconditioner::AMG_Preconditioner(const AMG_Preconditioner& other) : Preconditioner(PreconditionerType::AMG, other.pA), 
	amg_levels(other.amg_levels), niters_smooth(other.niters_smooth), w_cycle(other.w_cycle),
	Nm(other.Nm), Cm(other.Cm), Wm(other.Wm)
{
	Am.resize(amg_levels);
	for(int i = 0; i < amg_levels; i++)
		Am[i] = new SparseMatrix( *(other.Am[i]) );
}

AMG_Preconditioner& AMG_Preconditioner::operator=(const AMG_Preconditioner& other)
{
	mytype = PreconditionerType::AMG;
	pA = other.pA;
	amg_levels = other.amg_levels;
	niters_smooth = other.niters_smooth;
	w_cycle = other.w_cycle;
	Nm = other.Nm;
	Cm = other.Cm;
	Wm = other.Wm;
	for(int i = 0; i < amg_levels; i++)
		if(Am[i] != nullptr)	delete Am[i];
	Am.resize(amg_levels);
	for(int i = 0; i < amg_levels; i++)
		Am[i] = new SparseMatrix( *(other.Am[i]) );
	return *this;
}

AMG_Preconditioner::~AMG_Preconditioner() 
{
	for(int i = 0; i < amg_levels; i++)
		if(Am[i] != nullptr)	delete Am[i];
}

bool AMG_Preconditioner::SetupPreconditioner()
{
	return SetupAMG(0);
}

bool AMG_Preconditioner::PreconditionedSolve(const std::vector<double>& rhs, std::vector<double>& x)
{
	std::fill(x.begin(), x.end(), 0.0);
	if(w_cycle) 
		W_cycle(0, x, rhs);
	else 
		V_cycle(0, x, rhs);
	return true;
}

// prolongation from m to m-1, m > 0
void AMG_Preconditioner::prolongation(int m_from, const std::vector<double>& e, std::vector<double>& v_out)
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
void AMG_Preconditioner::restriction(int m_from, const std::vector<double>& e, std::vector<double>& v_out)
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
		int i = C[k];
		invC[i] = k;
		v_out[k] = e[i];
		processed[i] = true;
	}
	for(int j = 0; j < nf; j++) if(!processed[j])
	{
		const sparse_row& r = W[j];
		for(int k = 0; k < r.row.size(); k++)
		{
			int col = r.row[k].first;
			v_out[invC[col]] += r.row[k].second*e[j];
		}
	}
}

void AMG_Preconditioner::smoothing(int m, std::vector<double>& x, const std::vector<double>& b, int nrelax)
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
void AMG_Preconditioner::construct_coarse_system(int m_from)
{
	const SparseMatrix* Amf = Am[m_from];	// fine system

	int m = m_from+1;
	int nc = Nm[m];					// coarse system size
	int nf = Am[m_from]->Size();		// fine system size

	// (1) A^{m+1} = I_m^{m+1} A^m
	// (2) A^{m+1} = A^{m+1} I_{m+1}^m
	
	std::cout << "Fine system level " << m_from << " size " << Amf->Size() << std::endl;
	std::cout << "Coarse system level " << m << " size " << nc << std::endl;


	SparseMatrix* Amc = Am[m];
	Amc = new SparseMatrix(nc);

	// need to construct inverse mapping from C[k] to k to fill coarse system
	std::map<int,int> invC;

	// (1) A^{m+1} = I_m^{m+1} A^m
	const std::vector<int>& C = Cm[m];
	assert(nc == C.size());
	interpolation_type& W = Wm[m];
	std::vector<bool> coarse(nf,false);
	// coarse indices 
	for(int k = 0; k < C.size(); k++)
	{
		int i = C[k];
		invC[i] = k;
		coarse[i] = true;
	}
	for(int k = 0; k < C.size(); k++)
		Amc->set_vector(k, (*Amf)[C[k]], 0, nf);
	// fine indices
	for(int j = 0; j < nf; j++) if(!coarse[j])
	{
		const sparse_row& r = W[j];
		for(int k = 0; k < r.row.size(); k++)
		{
			int col = r.row[k].first;
			Amc->add_vector(invC[col], r.row[k].second, (*Amf)[j], 0, nf);
		}
	}

	// (2) A^{m+1} = A^{m+1} I_{m+1}^m
	// temporary intermediate matrix 
	SparseMatrix* Atmp = new SparseMatrix(nc);	
	for(int Arow = 0; Arow < nc; Arow++)
	{
		const sparse_row& r = (*Amc)[Arow];
		sparse_row& r2 = (*Atmp)[Arow];
		for(int k = 0; k < r.row.size(); k++)
		{
			int Acol = r.row[k].first;
			double Aelem = r.row[k].second;
			if(coarse[Acol])
				Atmp->add_element(Arow, invC[Acol], Aelem);
			else
			{
				const sparse_row& rW = W[Acol];
				for(int kk = 0; kk < rW.row.size(); kk++)
				{
					int col = rW.row[kk].first;
					if(fabs(rW.row[kk].second) > 1e-10)
						Atmp->add_element(Arow, invC[col], rW.row[kk].second * Aelem);
				}
			}
		}
		r2.remove_zeros();
	}

	std::swap(Amc, Atmp);
	delete Atmp;

	Am[m] = Amc;
}


void AMG_Preconditioner::construct_coarsest_inverse()
{
	return ILU0_Preconditioner::construct_inverse(Am[amg_levels-1]);
}

void AMG_Preconditioner::solve_coarsest(std::vector<double>& x, const std::vector<double>& b)
{
	return LU_in_place_solve(Am[amg_levels-1], b, x);
}

void find_strongly_connected(const SparseMatrix* pA, std::vector<std::vector<int> >& sc, double theta = 0.25)
{
	int n = pA->Size(), row, col;
	double val;
	std::vector<double> maxmod(n, 0.0);
	sc.resize(n);

	for(int row = 0; row < n; row++)
	{
		const sparse_row& r = (*pA)[row];
		for(int j = 0; j < r.row.size(); j++)
		{
			col = r.row[j].first;
			val = fabs(r.row[j].second);
			if(row != col && maxmod[row] < val)	maxmod[row] = val;
		}
	}

	for(int row = 0; row < n; row++)
	{
		const sparse_row& r = (*pA)[row];
		std::vector<int>& sci = sc[row];
		double rowmax = maxmod[row];
		for(int j = 0; j < r.row.size(); j++)
		{
			col = r.row[j].first;
			val = fabs(r.row[j].second);
			if(row != col && val > theta * rowmax)	sci.push_back(col);
		}
		// sort to use binary search over strong connections later
		if(!sci.empty())	std::sort(sci.begin(), sci.end());
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

void AMG_Preconditioner::construct_CF_partition(int m_from)
{
	bool debug = false;

	int m_to = m_from + 1;
	int m = m_to;

	const SparseMatrix* pmat = Am[m_from];
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
				q.ChangeKey(j, max_asc_amount);
				if( q.Pop() != j )	std::cerr << "Failed to remove " << j << " from queue" << std::endl, std::cin.ignore();
				const std::vector<int>& sc_j = strong_connections[j];
				for(int kk = 0; kk < sc_j.size(); kk++)
				{
					int l = sc_j[kk];
					if(q.Contains(l))	
						q.ChangeKey(l, q.GetKey(l)+1);
				}
			}
		}
		const std::vector<int>& sc_i = strong_connections[i];
		for(int k = 0; k < sc_i.size(); k++)
		{
			int j = sc_i[k];
			if(q.Contains(j))	
				q.ChangeKey(j, q.GetKey(j)-1);
		}
	}

	// extend constructed coarse partition	
	std::sort(C.begin(), C.end());
	std::sort(F.begin(), F.end());
	std::vector<bool> F_available(F.size(), true);
	int F_amount = F.size();	// amount of fine mesh elements

	if(debug)
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
		std::vector<int>::iterator Cb = C.begin(), Ce = C.end();
		// initialize connections
		for(int kk = 0; kk < sc_i.size(); kk++)
		{
			int l = sc_i[kk];
			if(std::binary_search(Cb, Ce, l))
				Im.push_back(l);
			else
				Ds.push_back(l);
		}
		std::sort(Im.begin(), Im.end());

		for(int kk = 0; kk < Ds.size(); kk++)
		{
			int j = Ds[kk];
			const std::vector<int>& sc_j = adj_strong_connections[j];//strong_connections[j];
			// intersect sc_j with Im \cup hat_Cm
			std::vector<int> v_intersection;
			std::set_intersection(sc_j.begin(), sc_j.end(), Im.begin(), Im.end(), std::back_inserter(v_intersection));
			if(v_intersection.empty() && !hat_Cm.empty())
				std::set_intersection(sc_j.begin(), sc_j.end(), hat_Cm.begin(), hat_Cm.end(), std::back_inserter(v_intersection));
			if(v_intersection.empty())
				hat_Cm.push_back(j);
		}
		std::vector<int>::iterator jt = Cb;
		if(hat_Cm.size() > 1 && !std::binary_search(Cb, Ce, i))
		{
			while(jt != Ce && i >= *jt)	++jt;
			C.insert(jt, i);
			F_available[k] = false;
			F_amount--;
		}
		else if(!hat_Cm.empty() && !std::binary_search(Cb, Ce, hat_Cm[0]))
		{
			int l = hat_Cm[0];
			while(jt != Ce && l >= *jt)	++jt;
			C.insert(jt, l);

			int kk = 0;
			while(kk < F.size() && F[kk] < l)	++kk;
			if(kk < F.size() && l == F[kk])
			{
				F_available[kk] = false;
				F_amount--;
			}
		}
	}

	if(debug)
	{
		std::cout << "after extension coarse " << C.size() << " fine " << F_amount << std::endl;
		for(int i = 0; i < C.size(); i++)
		{
			std::cout << C[i] << " ";
		}
		std::cout << std::endl;
		int nfine = 0;
		for(int i = 0; i < F.size(); i++) if(F_available[i])	nfine++;
		std::cout << "after extension fine " << nfine << std::endl;
		for(int i = 0; i < F.size(); i++) if(F_available[i])
		{
			std::cout << F[i] << " ";
		}
		std::cout << std::endl;
		std::cin.ignore();
	}

	Nm[m] = C.size();

	// repartition F so that disabled elements appear at the end
	int nF = F.size();
	int index1 = 0, index2 = nF-1;
	while(index1 < index2)
	{
		while(index1 < index2 &&  F_available[index1])	index1++;
		while(index1 < index2 && !F_available[index2])	index2--;
		if(index1 < index2)
		{
			std::swap(F[index1], F[index2]);
			F_available[index1] = true;
			F_available[index2] = false;
		}
		else if(index1 == nF-1 && F_available[index1])	index1++;
	}
	// remove disabled elements and sort 
	nF = index1;
	F.resize(nF);
	std::sort(F.begin(), F.end());

	assert(C.size() + F.size() == n);

	// determine interpolation weights
	interpolation_type& W = Wm[m];

	std::vector<int>::iterator Cb = C.begin(), Ce = C.end();
	int col;
	double val, di, sj, dk, aij_ajk;
	for(int k = 0; k < F.size(); k++)
	{
		int i = F[k];
		std::vector<int> Im_col;
		std::vector<double> Im_coef;
		std::vector<int> Ds;	// skipped strong connections
		const std::vector<int>& sc_i = strong_connections[i];
		// initialize connections
		// strong_connections is sorted, that's why Im and Ds are also sorted 
		const sparse_row& ri = (*pmat)[i];
		for(int kk = 0; kk < sc_i.size(); kk++)
		{
			int l = sc_i[kk];
			if(std::binary_search(Cb, Ce, l))
			{
				// find element with column l
				for(int kkk = 0; kkk < ri.row.size(); kkk++)
					if(ri.row[kkk].first == l)
					{
						Im_col.push_back(l);
						Im_coef.push_back(ri.row[kkk].second);
						break;
					}
			}
			else
				Ds.push_back(l);
		}
		
		if(debug)
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
			if(col == i || !std::binary_search(sc_i.begin(), sc_i.end(), col))
				di += val;
		}
		for(int kk = 0; kk < Ds.size(); kk++)
		{
			int j = Ds[kk];
			const std::vector<int>& sc_j = strong_connections[j];
			const sparse_row& rj = (*pmat)[j];
			sj = aij_ajk = 0.0;
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
					columns_values_Im_cap_Sj[col] += val*ajk;
				}
			}
			if(sj != 0.0)
			{
				std::vector<int>::const_iterator it_beg = Im_col.begin(), it_end = Im_col.end();
				for(std::map<int,double>::const_iterator itt = columns_values_Im_cap_Sj.begin(); itt != columns_values_Im_cap_Sj.end(); itt++)
				{
					std::vector<int>::const_iterator it = std::find(it_beg,it_end,itt->first);
					assert(it != it_end);
					int ind = it-it_beg;
					Im_coef[ind] += itt->second / sj;
				}
			}
		}
		if(!Im_col.empty())
		{
			sparse_row& rW = W[i];
			for(int kk = 0; kk < Im_col.size(); kk++)
				rW.add_element(std::make_pair(Im_col[kk], -Im_coef[kk] / di));
		}
		else
		{
			//?
		}
	}
}


void AMG_Preconditioner::V_cycle(int m, std::vector<double>& x, const std::vector<double>& b)
{
	const SparseMatrix* pA = Am[m];
	bool stop = m >= amg_levels-2;
	// pre-smoothing
	smoothing(m,x,b, niters_smooth);
	// compute residual r = b - Ax
	int nc = Am[m+1]->Size();
	int n = pA->Size();
	std::vector<double> resid(n), resid_restricted(nc), x0(nc, 0.0);
	for(int i = 0; i < n; i++)
	{
		const sparse_row& r = (*pA)[i];
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
		smoothing(m+1,x0,resid_restricted, niters_smooth);
#endif
	}	
	prolongation(m+1,x0,resid);
	for(int i = 0; i < n; i++)
		x[i] += resid[i];
	// post-smoothing
	smoothing(m,x,b, niters_smooth);
}

void AMG_Preconditioner::W_cycle(int m, std::vector<double>& x, const std::vector<double>& b)
{
	const SparseMatrix* pA = Am[m];
	bool stop = m >= amg_levels-2;
	// pre-smoothing
	smoothing(m,x,b, niters_smooth);
	// compute residual r = b - Ax
	int nc = Am[m+1]->Size();
	int n = pA->Size();
	std::vector<double> resid(n), resid_restricted(nc), x0(nc, 0.0);
	for(int i = 0; i < n; i++)
	{
		const sparse_row& r = (*pA)[i];
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
		const sparse_row& r = (*pA)[i];
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

bool AMG_Preconditioner::SetupAMG(int m)
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




#if 0
void AMG_Preconditioner::construct_CF_partition(int m_from)
{
	bool debug = false;

	int m_to = m_from + 1;
	int m = m_to;

	const SparseMatrix* pmat = Am[m_from];
	int n = pmat->Size();
	PriorityQueue<int, std::greater<int> > q(n);
	std::vector< std::vector<int> > strong_connections(n);		// strong connections for each grid element
	std::vector< std::vector<int> > adj_strong_connections(n);	// adjacent strong connections for each element

	std::vector<int>& C = Cm[m];	// coarse	partition on level m
	//std::vector<int> F;				// fine 	partition on level m

	std::set<int> C2;
	std::set<int> F;

	find_strongly_connected(pmat, strong_connections);
	find_adjacent_strongly_connected(strong_connections, adj_strong_connections);
	for(int i = 0; i < adj_strong_connections.size(); i++)
		q.Push(i, adj_strong_connections[i].size());

	int max_asc_amount = n;	// maximum amount of adjacent strong connections cannot be more than matrix size
	while(!q.Empty())
	{
		int i = q.Peek();
		//C.push_back(i);
		q.Pop();
		const std::vector<int>& asc_i = adj_strong_connections[i];
		for(int k = 0; k < asc_i.size(); k++)
		{
			int j = asc_i[k];
			if(q.Contains(j))
			{
				F.insert(j);
				// remove j from q
				q.ChangeKey(j, max_asc_amount);
				if( q.Pop() != j )	std::cerr << "Failed to remove " << j << " from queue" << std::endl, std::cin.ignore();
				const std::vector<int>& sc_j = strong_connections[j];
				for(int kk = 0; kk < sc_j.size(); kk++)
				{
					int l = sc_j[kk];
					if(q.Contains(l))	
						q.ChangeKey(l, q.GetKey(l)+1);
				}
			}
		}
		const std::vector<int>& sc_i = strong_connections[i];
		for(int k = 0; k < sc_i.size(); k++)
		{
			int j = sc_i[k];
			if(q.Contains(j))	
				q.ChangeKey(j, q.GetKey(j)-1);
		}
	}

	std::sort(C.begin(), C.end());
	//std::sort(F.begin(), F.end());

	if(debug)
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

	// extend constructed coarse partition	
	//std::vector<bool> F_available(F.size(), true);
	//int F_amount = F.size();	// amount of fine mesh elements

	for(int k = 0; k < F.size(); k++)
	{
		if(!F_available[k])	continue;	// skip disabled elements 
		int i = F[k];
		std::set<int> hat_Cm;
		std::set<int> Im;	// interpolation connections
		std::vector<int> Ds;	// skipped strong connections
		const std::vector<int>& sc_i = strong_connections[i];
		// initialize connections
		for(int kk = 0; kk < sc_i.size(); kk++)
		{
			int l = sc_i[kk];
			if(std::binary_search(C.begin(), C.end(), l))
				Im.insert(l);
			else
				Ds.push_back(l);
		}
		std::sort(Im.begin(), Im.end());

		for(int kk = 0; kk < Ds.size(); kk++)
		{
			int j = Ds[kk];
			const std::vector<int>& sc_j = adj_strong_connections[j];//strong_connections[j];
			// intersect sc_j with Im \cup hat_Cm
			std::vector<int> v_intersection;
			std::set_intersection(sc_j.begin(), sc_j.end(), Im.begin(), Im.end(), std::back_inserter(v_intersection));
			if(!hat_Cm.empty() && v_intersection.empty())
				std::set_intersection(sc_j.begin(), sc_j.end(), hat_Cm.begin(), hat_Cm.end(), std::back_inserter(v_intersection));
			if(v_intersection.empty())
			{
				//std::vector<int>::iterator jt = Im.begin();
				//while(jt != Im.end() && j >= *jt)	++jt;
				//Im.insert(jt, j);
				hat_Cm.insert(j);//hat_Cm.push_back(j);
			}
		}
		std::vector<int>::iterator Cb = C.begin(), Ce = C.end(), jt = Cb;
		if(hat_Cm.size() > 1 && !std::binary_search(Cb, Ce, i))
		{
			while(jt != C.end() && i >= *jt)	++jt;
			C.insert(jt, i);
			// if(F_available[k])
			// {
			// 	F_available[k] = false;
			// 	F_amount--;
			// }
		}
		else if(!hat_Cm.empty() && !std::binary_search(Cb, Ce, *(hat_Cm.begin())/*hat_Cm[0]*/))
		{
			int l = *(hat_Cm.begin());//hat_Cm[0];
			while(jt != Ce && l >= *jt)	++jt;
			C.insert(jt, l);

			// int kk = 0;
			// while(kk < F.size() && F[kk] < l)	++kk;
			// if(kk < F.size() && l == F[kk])
			// {
			// 	F_available[kk] = false;
			// 	F_amount--;
			// }
		}
	}

	if(debug)
	{
		std::cout << "after extension coarse " << C.size() << " fine " << F_amount << std::endl;
		for(int i = 0; i < C.size(); i++)
		{
			std::cout << C[i] << " ";
		}
		std::cout << std::endl;
		int nfine = 0;
		for(int i = 0; i < F.size(); i++) if(F_available[i])	nfine++;
		std::cout << "after extension fine " << nfine << std::endl;
		for(int i = 0; i < F.size(); i++) if(F_available[i])
		{
			std::cout << F[i] << " ";
		}
		std::cout << std::endl;
		std::cin.ignore();
	}

	Nm[m] = C.size();

	// repartition F so that disabled elements appear at the end
	int nF = F.size();
	int index1 = 0, index2 = nF-1;
	while(index1 < index2)
	{
		while(index1 < index2 &&  F_available[index1])	index1++;
		while(index1 < index2 && !F_available[index2])	index2--;
		if(index1 < index2)
		{
			std::swap(F[index1], F[index2]);
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
		std::vector<int> Im_col;
		std::vector<double> Im_coef;
		std::vector<int> Ds;	// skipped strong connections
		const std::vector<int>& sc_i = strong_connections[i];
		// initialize connections
		// strong_connections is sorted, that's why Im and Ds are also sorted 
		const sparse_row& ri = (*pmat)[i];
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
					std::swap(Im_col[kk], Im_col[kkk]);
					std::swap(Im_coef[kk], Im_coef[kkk]);
				}
			}
		
		if(debug)
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
			if(col == i || !std::binary_search(sc_i.begin(), sc_i.end(), col))
				di += val;
		}
		for(int kk = 0; kk < Ds.size(); kk++)
		{
			int j = Ds[kk];
			const std::vector<int>& sc_j = strong_connections[j];
			const sparse_row& rj = (*pmat)[j];
			sj = 0.0;
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
					columns_values_Im_cap_Sj[col] += val*ajk;
				}
			}
			if(sj != 0.0)
			{
				std::vector<int>::const_iterator it_beg = Im_col.begin(), it_end = Im_col.end();
				for(std::map<int,double>::const_iterator itt = columns_values_Im_cap_Sj.begin(); itt != columns_values_Im_cap_Sj.end(); itt++)
				{
					std::vector<int>::const_iterator it = std::find(it_beg,it_end,itt->first);
					assert(it != it_end);
					int ind = it-it_beg;
					Im_coef[ind] += itt->second / sj;
				}
			}
		}
		if(!Im_col.empty())
		{
			sparse_row& rW = W[i];
			for(int kk = 0; kk < Im_col.size(); kk++)
				rW.add_element(std::make_pair(Im_col[kk], -Im_coef[kk] / di));
		}
		else
		{
			//?
		}
	}
}
#endif