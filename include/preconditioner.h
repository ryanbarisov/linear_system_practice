#ifndef _PRECONDITIONER_H
#define _PRECONDITIONER_H

#include <matrix.h>
#include "solver_parameters.hpp"
#include "priority_queue.hpp"


enum class PreconditionerType
{
	NONE,
	ILU0,
	ILUC,
	AMG,
	JACOBI,
	GAUSS_SEIDEL
};

class Preconditioner
{
protected:
	PreconditionerType mytype;
	const SparseMatrix* pA;
public:
	Preconditioner(PreconditionerType _mytype, const SparseMatrix* _pA) : mytype(_mytype), pA(_pA) {}
	PreconditionerType GetType() const {return mytype;}

	virtual ~Preconditioner() {}
	
	virtual bool SetupPreconditioner() = 0;	
	virtual bool PreconditionedSolve(const std::vector<double>& rhs, std::vector<double>& x) = 0;
};

class ILU0_Preconditioner : public Preconditioner
{
private:
	SparseMatrix * LU;
public:
	ILU0_Preconditioner(const SparseMatrix* _pA, SolverParameters parameters);
	~ILU0_Preconditioner();

	bool SetupPreconditioner() override;
	bool PreconditionedSolve(const std::vector<double>& rhs, std::vector<double>& x) override;

	static void construct_inverse(SparseMatrix* pA);
};

class ILUC_Preconditioner : public Preconditioner
{
private:
	SparseMatrix * L;
	SparseMatrix * U;

	double tau;
	int lfil;
public:
	ILUC_Preconditioner(const SparseMatrix* _pA, SolverParameters parameters);
	~ILUC_Preconditioner();

	bool SetupPreconditioner() override;
	bool PreconditionedSolve(const std::vector<double>& rhs, std::vector<double>& x) override;
};


#if !defined(RECURSIVE_AMG)

class AMG_Preconditioner : public Preconditioner
{
	typedef std::map<int,sparse_row> interpolation_type;
private:
	int amg_levels = 2;
	int niters_smooth = 2;
	int w_cycle = 0;

	std::vector<int> Nm;						// coarse system dimensions for 1 <= m <= q = amg_levels-1
	std::vector< std::vector<int> > Cm;			// index set for coarse grids, 1 <= m <= q = amg_levels-1
	std::vector< SparseMatrix* > Am;				// coarse grid matrices for 0 <= m <= amg_levels-1
	std::vector< interpolation_type > Wm;		// interpolation weights for 1 <= m <= amg_levels-1

	// prolongation from m to m-1, m > 0
	void prolongation(int m_from, const std::vector<double>& e, std::vector<double>& v_out);
	// restriction from m to m+1
	void restriction(int m_from, const std::vector<double>& e, std::vector<double>& v_out);
	void smoothing(int m, std::vector<double>& x, const std::vector<double>& b, int nrelax = 1);

	void construct_coarsest_inverse();
	void solve_coarsest(std::vector<double>& x, const std::vector<double>& b);

	// construct system on step m+1 with given system on step m and interpolator matrices
	// A^{m+1} = I_m^{m+1} A^m I_{m+1}^m
	void construct_coarse_system(int m_from);
	void construct_CF_partition(int m_from);

	bool SetupAMG(int lvl = 0);

	void V_cycle(int m, std::vector<double>& x, const std::vector<double>& b);
	void W_cycle(int m, std::vector<double>& x, const std::vector<double>& b);

public:
	AMG_Preconditioner(const SparseMatrix* _pA, SolverParameters parameters);
	~AMG_Preconditioner();

	bool SetupPreconditioner() override;
	bool PreconditionedSolve(const std::vector<double>& rhs, std::vector<double>& x) override;
};

#else

class AMG_Preconditioner : public Preconditioner
{
	Preconditioner* Next = nullptr;
public:
	AMG_Preconditioner(const SparseMatrix* _pA, SolverParameters parameters);
	AMG_Preconditioner(const AMG_Preconditioner* prev, int amg_levels, int depth);
	~AMG_Preconditioner();

	bool SetupPreconditioner() override;
	bool PreconditionedSolve(const std::vector<double>& rhs, std::vector<double>& x) override;
};

#endif

class JACOBI_Preconditioner : public Preconditioner
{
public:
	JACOBI_Preconditioner(const SparseMatrix* _pA, SolverParameters parameters) : Preconditioner(PreconditionerType::JACOBI, _pA) {}

	bool SetupPreconditioner() {return true;}
	bool PreconditionedSolve(const std::vector<double>& rhs, std::vector<double>& x) override;
};

class GS_Preconditioner : public Preconditioner
{
public:
	GS_Preconditioner(const SparseMatrix* _pA, SolverParameters parameters) : Preconditioner(PreconditionerType::GAUSS_SEIDEL, _pA) {}

	bool SetupPreconditioner() {return true;}
	bool PreconditionedSolve(const std::vector<double>& rhs, std::vector<double>& x) override;
};

////////////////////////////////////////////////////////////////////

static Preconditioner * CreatePreconditioner(PreconditionerType ptype, const SparseMatrix* pA, const SolverParameters& parameters)
{
	if(ptype == PreconditionerType::ILU0)
	{
		return new ILU0_Preconditioner(pA, parameters);
	}
	if(ptype == PreconditionerType::ILUC)
	{
		return new ILUC_Preconditioner(pA, parameters);
	}
	else if(ptype == PreconditionerType::AMG)
	{
		return new AMG_Preconditioner(pA, parameters);
	}
	else if(ptype == PreconditionerType::JACOBI)
	{
		return new JACOBI_Preconditioner(pA, parameters);
	}
	else if(ptype == PreconditionerType::GAUSS_SEIDEL)
	{
		return new GS_Preconditioner(pA, parameters);
	}
	else
	{
		std::cerr << "Sorry, preconditioner was not implemented " << std::endl;
		return NULL;
	}
}


#endif