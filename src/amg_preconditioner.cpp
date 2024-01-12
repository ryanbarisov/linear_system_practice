#include <preconditioner.h>
#include <solver_parameters.h>
#include <matrix.h>
#include <method.h>
#include <priority_queue.hpp>
#include <algorithm>
#include <set>

#define AMG_COARSE_SOLVE // solve coarsest grid instead of smoothing


static void find_strongly_connected(const CSRMatrix* pA, std::vector<std::vector<int> >& sc, double theta = 0.25)
{
	int n = pA->Size();
	std::vector<double> maxmod(n, 0.0);
	sc.resize(n);

	for(int i = 0; i < n; ++i)
	{
		double m = 0.0;
		for(int j = pA->GetIA(i); j < pA->GetIA(i+1); ++j)
			if(pA->GetJA(j) != i && fabs(pA->GetA(j)) > m)
				m = fabs(pA->GetA(j));
		maxmod[i] = m;
	}

	for(int i = 0; i < n; ++i)
	{
		double m = maxmod[i];
		std::vector<int>& sci = sc[i];
		for(int j = pA->GetIA(i); j < pA->GetIA(i+1); ++j)
			if(pA->GetJA(j) != i && fabs(pA->GetA(j)) > theta*m)
				sci.push_back(pA->GetJA(j));
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

// prolongation from m to m-1, m > 0
void AMG_Preconditioner::prolongation(int m, const std::vector<double>& e, std::vector<double>& v_out) const
{
	int n = Am[m-1]->Size();
	v_out.resize(n);
	std::fill(v_out.begin(), v_out.end(), 0.0);
	const std::vector<int>& C = Cm[m-1];
	const interpolation_type& W = Wm[m-1];
	assert(e.size() == C.size());
	for(int k = 0; k < C.size(); k++)
	{
		int i = C[k];
		v_out[i] = e[k];
		if(W.find(i) != W.end())
		{
			const RowAccumulator& r = W.at(i);
			for(int j = 0; j < r.w.size(); ++j)
			{
				assert(r.jw[j] != i);
				v_out[r.jw[j]] += r.w[j]*e[k];
			}
		}
	}
}
// restriction from m to m+1
void AMG_Preconditioner::restriction(int m, const std::vector<double>& e, std::vector<double>& v_out) const
{
	int nf = e.size();
	int nc = Am[m+1]->Size();
	v_out.resize(nc);
	const std::vector<int>& C = Cm[m];
	const interpolation_type& W = Wm[m];

	for(int k = 0; k < C.size(); k++)
	{
		int i = C[k];
		v_out[k] = e[i];
		if(W.find(i) != W.end())
		{
			const RowAccumulator& r = W.at(i);
			for(int j = 0; j < r.w.size(); ++j)
			{
				assert(r.jw[j] != i);
				v_out[k] += r.w[j]*e[r.jw[j]];
			}
		}
	}
}

void AMG_Preconditioner::smoothing(int m, std::vector<double>& x, const std::vector<double>& b, int nrelax)
{
	int method = 1;	// 0 - Jacobi, 1 - Gauss-Seidel, 2 - ?
	for(int k = 0; k < nrelax; k++)
	{
		if(method == 0)
			jacobi_precondition(Am[m], x, b);
		else if(method == 1)
			gs_precondition(Am[m], x, b);
	}
}

// construct system on step m+1 with given system on step m and interpolator matrices
// A^{m+1} = I_m^{m+1} A^m I_{m+1}^m
void AMG_Preconditioner::construct_coarse_system(int m)
{
	const CSRMatrix* pAf = Am[m-1];	// fine system

	int nf = pAf->Size();// fine system size
	int nc = Nm[m];	// coarse system size
	int nfc= nf-nc;

	// (1) A^{m+1} = I_m^{m+1} A^m
	// (2) A^{m+1} = A^{m+1} I_{m+1}^m
	
	std::cout << "Fine system level " << m-1 << " size " << nf << std::endl;
	std::cout << "Coarse system level " << m << " size " << nc << std::endl;

	const std::vector<int>& C = Cm[m-1];
	assert(nc == C.size());
	const interpolation_type& W = Wm[m-1];
	// (1) A^{m+1} = I_m^{m+1} A^m
	RowAccumulator w(nf);
	std::vector<int> ib(nc+1), jb;
	std::vector<double> b;
	for(int k = 0; k < nc; ++k)
	{
		int i = C[k];
		w.SetRowFrom(pAf, i);
		//assert(W.find(i) != W.end());
		if(W.find(i) != W.end())
		{
			const RowAccumulator& rW = W.at(i);
			for(int j = 0; j < rW.jw.size(); ++j)
				w.SparseAdd(pAf, rW.jw[j], rW.w[j]);
		}
		ib[k+1] = ib[k] + w.Size();
		for(int k = 0; k < w.Size(); ++k)
		{
			jb.push_back(w.jw[k]);
			b.push_back(w.w[k]);
		}
		w.Clear();
	}
	CSRMatrix* pAc = new CSRMatrix(b,ib,jb);
	pAc->RemoveZeros();
	pAc->FlipStorageFormat(nf); // convert to CSC
	ib.resize(nc+1); ib[0] = 0;
	b.clear();jb.clear();
	// (2) A^{m+1} = A^{m+1} I_{m+1}^m
	for(int k = 0; k < nc; ++k)
	{
		int i = C[k];
		w.SetRowFrom(pAc, i);
		//assert(W.find(i) != W.end());
		if(W.find(i) != W.end())
		{
			const RowAccumulator& rW = W.at(i);
			for(int j = 0; j < rW.jw.size(); ++j)
				w.SparseAdd(pAc, rW.jw[j], rW.w[j]);
		}
		ib[k+1] = ib[k] + w.Size();
		for(int k = 0; k < w.Size(); ++k)
		{
			jb.push_back(w.jw[k]);
			b.push_back(w.w[k]);
		}
		w.Clear();
	}
	delete pAc;
	Am[m] = new CSRMatrix(b,ib,jb, false);
	Am[m]->RemoveZeros();
	Am[m]->FlipStorageFormat(nc); // back to CSR
	//std::string s = "amg_level_"+std::to_string(m)+".csr";
	//Am[m]->Save(s.c_str(), MatrixFormat::CSR);
}


void AMG_Preconditioner::construct_coarsest_inverse()
{
	return ILU0_Preconditioner::construct_inverse(Am[amg_levels-1]);
}

void AMG_Preconditioner::solve_coarsest(std::vector<double>& x, const std::vector<double>& b)
{
	return LU_in_place_solve(Am[amg_levels-1], b, x);
}

void AMG_Preconditioner::construct_CF_partition(int m)
{
	bool debug = false;

	const CSRMatrix* pA = Am[m-1]; // fine system matrix
	int n = pA->Size();
	PriorityQueue<int, std::greater<int> > q(n);
	std::vector< std::vector<int> > strong_connections(n);		// strong connections for each grid element
	std::vector< std::vector<int> > adj_strong_connections(n);	// adjacent strong connections for each element

	std::vector<int>& C = Cm[m-1];	// coarse partition on level m
	std::vector<int> F;		// fine partition on level m

	find_strongly_connected(pA, strong_connections);
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
		//std::sort(Im.begin(), Im.end());

		for(int kk = 0; kk < Ds.size(); kk++)
		{
			int j = Ds[kk];
			const std::vector<int>& sc_j = adj_strong_connections[j];
			// intersect sc_j with Im V hat_Cm
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
	interpolation_type& W = Wm[m-1];

	std::vector<int>::iterator Cb = C.begin(), Ce = C.end();
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
		for(int kk = 0; kk < sc_i.size(); kk++)
		{
			int l = sc_i[kk];
			if(std::binary_search(Cb, Ce, l))
			{
				// find element with column l
				for(int j = pA->GetIA(i); j < pA->GetIA(i+1); ++j)
					if(pA->GetJA(j) == l)
					{
						Im_col.push_back(l);
						Im_coef.push_back(pA->GetA(j));
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
		for(int j = pA->GetIA(i); j < pA->GetIA(i+1); ++j)
		{
			col = pA->GetJA(j);
			if(col == i || !std::binary_search(sc_i.begin(), sc_i.end(), col))
				di += pA->GetA(j);
		}
		for(int kk = 0; kk < Ds.size(); kk++)
		{
			int j = Ds[kk];
			const std::vector<int>& sc_j = strong_connections[j];
			sj = 0.0;
			std::map<int,double> columns_values_Im_cap_Sj;
			std::vector<int> intersection;	// work vector to store sorted array intersections
			std::set_intersection(sc_j.begin(),sc_j.end(),Im_col.begin(),Im_col.end(),std::back_inserter(intersection));
			for(int j1 = pA->GetIA(j); j1 < pA->GetIA(j+1); ++j1)
			{
				col = pA->GetJA(j1);
				if(std::binary_search(intersection.begin(),intersection.end(),col))
				{
					double ajk = pA->GetA(j1);
					sj += ajk;
					columns_values_Im_cap_Sj[col] += pA->Get(i,j)*ajk;
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
			/*
			RowAccumulator& rW = W[i];
			rW.Resize(n);
			for(int kk = 0; kk < Im_col.size(); kk++)
				rW.Add(Im_col[kk], -Im_coef[kk] / di);
			*/
			int nfc = nF;
			for(int kk = 0; kk < Im_col.size(); kk++)
				if(W[Im_col[kk]].Size() != n)
					W[Im_col[kk]].Resize(n);
			for(int kk = 0; kk < Im_col.size(); kk++)
				W[Im_col[kk]].Push(i, -Im_coef[kk] / di);
		}
		else
		{
			//?
		}
	}
}


AMG_Preconditioner::AMG_Preconditioner(const CSRMatrix* pA, const SolverParameters& params)
	: Preconditioner(pA, params)
	{
		amg_levels = std::max(2, params.GetIntegerParameter("amg_levels").second);
		niters_smooth = std::max(0, params.GetIntegerParameter("niters_smooth").second);
		w_cycle = std::max(0, params.GetIntegerParameter("w_cycle").second);
		Nm.resize(amg_levels);
		Cm.resize(amg_levels-1);
		Am.resize(amg_levels);
		Wm.resize(amg_levels-1);
		Am[0] = new CSRMatrix(pA->GetA(), pA->GetIA(), pA->GetJA());
		Nm[0] = Am[0]->Size();
		//for(int i = 0; i < Nm[0]; i++)	Cm[0].push_back(i);	// is it needed?
	}

AMG_Preconditioner::~AMG_Preconditioner()
{
	for(int i = 0; i < amg_levels; i++)
		if(Am[i] != nullptr)	delete Am[i];
}

bool AMG_Preconditioner::SetupPreconditioner()
{
	return SetupAMG();
}

bool AMG_Preconditioner::SetupAMG(int m)
{
	if(m+1 < amg_levels)
	{
		construct_CF_partition(m+1);
		construct_coarse_system(m+1);
		SetupAMG(m+1);
	}
#if defined(AMG_COARSE_SOLVE)
	else
		construct_coarsest_inverse();
#endif
	return true;
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


void AMG_Preconditioner::V_cycle(int m, std::vector<double>& x, const std::vector<double>& b)
{
	const CSRMatrix* pA = Am[m];
	bool stop = m >= amg_levels-2;
	// pre-smoothing
	smoothing(m,x,b, niters_smooth);
	
	int nc = Am[m+1]->Size();
	int n = pA->Size();
	std::vector<double> resid = b, resid_restricted(nc), x0(nc, 0.0);
	pA->Multiply(-1.0,x,1.0,resid); // r = b - Ax
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
	for(int i = 0; i < n; i++) x[i] += resid[i];
	// post-smoothing
	smoothing(m,x,b, niters_smooth);
}

void AMG_Preconditioner::W_cycle(int m, std::vector<double>& x, const std::vector<double>& b)
{
	const CSRMatrix* pA = Am[m];
	bool stop = m >= amg_levels-2;
	// pre-smoothing
	smoothing(m,x,b, niters_smooth);
	// compute residual r = b - Ax
	int nc = Am[m+1]->Size();
	int n = pA->Size();
	std::vector<double> resid = b, resid_restricted(nc), x0(nc, 0.0);
	pA->Multiply(-1.0, x, 1.0,resid);
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
	resid = b;
	pA->Multiply(-1.0, x, 1.0,resid);
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

