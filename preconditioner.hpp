#ifndef _PRECONDITIONER_H
#define _PRECONDITIONER_H

#include "matrix.hpp"

#include <map>
#include <utility>

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