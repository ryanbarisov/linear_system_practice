#include <algorithm>
#include <preconditioner.h>
#include <solver_parameters.h>
#include <matrix.h>
#include <method.h>

double Timer();

ILUC_Preconditioner::ILUC_Preconditioner(const CSRMatrix* pA, const SolverParameters& params)
	: Preconditioner(pA, params), L(nullptr), U(nullptr)
	{
		tau = params.GetRealParameter("drop_tolerance").second;
		lfil = params.GetIntegerParameter("level_of_fill").second;
	}

ILUC_Preconditioner::~ILUC_Preconditioner() 
{
	if(L != nullptr)	delete L;
	if(U != nullptr)	delete U;
}


//https://faculty.cc.gatech.edu/~echow/pubs/crout.pdf
bool ILUC_Preconditioner::SetupPreconditioner()
{
	int n = pA->Size();
	std::vector<int> Ufirst(n), Lfirst(n);
	std::vector<std::vector<int>> Ulist(n), Llist(n);
	std::vector<std::vector<int>> Alist(n);
	std::vector<std::vector<double>> Blist(n);
	for(int i = 0; i < n; ++i)
		for(int j = pA->GetIA(i); j < pA->GetIA(i+1); ++j)
			if(pA->GetJA(j) < i)
			{
				Alist[pA->GetJA(j)].push_back(i);
				Blist[pA->GetJA(j)].push_back(pA->GetA(j));
			}
	RowAccumulator z(n), w(n);
	std::vector<int> ila(n+1), jla, iua(n+1), jua;
	std::vector<double> la, ua;
	jla.reserve(n*(2*lfil+1));
	jua.reserve(n*(2*lfil+1));
	la.reserve(n*(2*lfil+1));
	ua.reserve(n*(2*lfil+1));
	double t_sadd = 0.0, t_upd1 = 0.0, t_upd2 = 0.0, t_clr = 0.0, t0;
	for(int k = 0; k < n; ++k)
	{
		t0 = Timer();
		z.SetRowFrom(pA, k, k, n);
		double znorm = z.L2Norm();
		std::vector<int>& Ll = Llist[k], &Ul = Ulist[k];
		for(int j = 0; j < Ll.size(); ++j)
		{
			int i = Ll[j];
			assert(jla[Lfirst[i]] == k);
			z.SparseAdd(ua,jua,-la[Lfirst[i]],Ufirst[i],iua[i+1]);
		}
		w.SetIntervalFrom(pA->Size(), Alist[k], Blist[k]);
		double wnorm = w.L2Norm();
		for(int j = 0; j < Ul.size(); ++j)
		{
			int i = Ul[j];
			assert(jua[Ufirst[i]] == k);
			w.SparseAdd(la,jla,-ua[Ufirst[i]],Lfirst[i],ila[i+1]);
		}
		t_sadd += Timer()-t0;
		z.Drop(k,tau*znorm,lfil);
		w.Drop(k,tau*wnorm,lfil);
		w.Remove(k);
		Ufirst[k] = iua[k];
		Lfirst[k] = ila[k];
		t0 = Timer();
		std::copy(z.jw.begin(), z.jw.end(), std::back_inserter(jua));
		std::copy(z.w.begin(), z.w.end(), std::back_inserter(ua));
		iua[k+1] = iua[k] + z.Size();
		for(int j = 0; j < Ul.size(); ++j)
		{
			int i = Ul[j];
			if(Ufirst[i]+1 < iua[i+1])
				Ulist[jua[++Ufirst[i]]].push_back(i);
		}
		Ulist[jua[++Ufirst[k]]].push_back(k);
		t_upd1 += Timer()-t0;

		t0 = Timer();
		ila[k+1] = ila[k] + w.Size()+1;
		jla.push_back(k);
		la.push_back(1.0);
		w.Scale(1.0/z.Get(k));
		std::copy(w.jw.begin(), w.jw.end(), std::back_inserter(jla));
		std::copy(w.w.begin(), w.w.end(), std::back_inserter(la));
		for(int j = 0; j < Ll.size(); ++j)
		{
			int i = Ll[j];
			if(Lfirst[i]+1 < ila[i+1])
				Llist[jla[++Lfirst[i]]].push_back(i);
		}
		Llist[jla[++Lfirst[k]]].push_back(k);
		t_upd2 += Timer()-t0;

		t0 = Timer();
		z.Clear();
		w.Clear();
		t_clr += Timer()-t0;
	}
	assert(ila[n] == la.size());
	L = new CSRMatrix(la,ila,jla);
	assert(iua[n] == ua.size());
	U = new CSRMatrix(ua,iua,jua);

	//std::cout << "Sparse add time: " << t_sadd << std::endl;
	//std::cout << "Update time: z " << t_upd1 << " w " << t_upd2 << std::endl;
	//std::cout << "Clear time: " << t_clr << std::endl;

	return true;
}


bool ILUC_Preconditioner::PreconditionedSolve(const std::vector<double>& rhs, std::vector<double>& x)
{
	LU_solve(*L,*U,rhs,x);
	return true;
}
