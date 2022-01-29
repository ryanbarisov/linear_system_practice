#include <preconditioner.h>
#include <solver_parameters.hpp>
#include <matrix.h>
#include <method.h>
#include <list>


ILUC_Preconditioner::ILUC_Preconditioner(const SparseMatrix* pA, SolverParameters parameters) 
	: Preconditioner(PreconditionerType::ILUC, pA), L(nullptr), U(nullptr) 
	{
		tau = parameters.GetRealParameter("drop_tolerance").second;
		lfil = parameters.GetIntegerParameter("level_of_fill").second;
	}

ILUC_Preconditioner::~ILUC_Preconditioner() 
{
	if(L != nullptr)	delete L;
	if(U != nullptr)	delete U;
}

bool ILUC_Preconditioner::SetupPreconditioner()
{
	int n = pA->Size();
	L = new SparseMatrix(n);
	U = new SparseMatrix(n);

	// prepare rows and columns of A into separate structures
	// sparse_row * Arows = new sparse_row[n];
	sparse_row * Acols = new sparse_row[n];
	entry e;
	for(int k = 0; k < pA->Size(); k++)
	{
		const sparse_row& r = (*pA)[k];
		// Arows[k] = r;
		for(int j = 0; j < r.row.size(); j++)
		{
			int col = r.row[j].first;
			Acols[col].add_element(k, r.row[j].second);
		}
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
		z.assign((*pA)[k],k,n);
		for(std::list<entry>::const_iterator it = Llist[k].begin(); it != Llist[k].end(); ++it)
			z.plus(-it->second, (*U)[it->first],k,n);
		w.assign(Acols[k],k+1,n);
		for(std::list<entry>::const_iterator it = Ulist[k].begin(); it != Ulist[k].end(); ++it)
			w.plus(-it->second, (*L)[it->first],k+1,n);
		
		double ukk = z.get_element(k);
		w *= 1.0/ukk;
		// TODO dropping rule for z and w

		// drop row z
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
		// drop column w
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
		//L->set_vector(k,w,k+1,n);

		for(sparse_type::iterator it = w.row.begin(); it != w.row.end(); ++it)
		{
			int row = it->first;
			double val = it->second;
			L->add_element(row,k,val);
		}

		L->add_element(k,k,1.0);

		// update indices
		for(int i = 0; i < k+1; i++)
		{
			const sparse_row& rU = (*U)[i];
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
			const sparse_row& rL = (*L)[i];
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
	// delete [] Arows;
	delete [] Acols;

	return true;
}

bool ILUC_Preconditioner::PreconditionedSolve(const std::vector<double>& rhs, std::vector<double>& x)
{
	LU_solve(*L,*U,rhs,x);
	return true;
}