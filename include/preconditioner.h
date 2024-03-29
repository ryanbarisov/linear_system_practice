#ifndef _PRECONDITIONER_H
#define _PRECONDITIONER_H

#include <matrix.h>
#include <solver_parameters.h>


enum class PreconditionerType
{
	NONE,
	ILU0,
	ILUC,
	AMG,
	JACOBI,
	GAUSS_SEIDEL,
	SSOR
};


PreconditionerType GetPreconditionerType(std::string name);

class Preconditioner
{
protected:
	PreconditionerType mytype;
	const CSRMatrix* pA;
public:
	Preconditioner(const CSRMatrix* _pA, const SolverParameters& params) : pA(_pA)
	{
		mytype = GetPreconditionerType(params.GetStringParameter("preconditioner").second);
	}
	PreconditionerType GetType() const {return mytype;}

	virtual ~Preconditioner() {}
	
	virtual bool SetupPreconditioner() = 0;	
	virtual bool PreconditionedSolve(const std::vector<double>& rhs, std::vector<double>& x) = 0;
};

class NOOP_Preconditioner : public Preconditioner
{
public:
	NOOP_Preconditioner(const CSRMatrix* _pA, const SolverParameters& params) : Preconditioner(_pA, params) {}
	bool SetupPreconditioner() {return true;}
	bool PreconditionedSolve(const std::vector<double>& rhs, std::vector<double>& x)
	{
		std::copy(rhs.begin(), rhs.end(), x.begin());
		return true;
	}
};

class ILU0_Preconditioner : public Preconditioner
{
private:
	CSRMatrix * LU;
public:
	ILU0_Preconditioner(const CSRMatrix* _pA, const SolverParameters& params);
	~ILU0_Preconditioner();

	bool SetupPreconditioner() override;
	bool PreconditionedSolve(const std::vector<double>& rhs, std::vector<double>& x) override;

	static void construct_inverse(CSRMatrix* pA);
};

class ILUC_Preconditioner : public Preconditioner
{
private:
	CSRMatrix * L, * U;

	double tau;
	int lfil;
public:
	ILUC_Preconditioner(const CSRMatrix* _pA, const SolverParameters& params);
	~ILUC_Preconditioner();

	bool SetupPreconditioner() override;
	bool PreconditionedSolve(const std::vector<double>& rhs, std::vector<double>& x) override;
};


#if !defined(RECURSIVE_AMG)

class AMG_Preconditioner : public Preconditioner
{
	typedef std::map<int,RowAccumulator> interpolation_type;
private:
	int amg_levels = 2;
	int niters_smooth = 2;
	int w_cycle = 0;

	std::vector<int> Nm;						// coarse system dimensions for 1 <= m <= q = amg_levels-1
	std::vector< std::vector<int> > Cm;			// index set for coarse grids, 1 <= m <= q = amg_levels-1
	std::vector< CSRMatrix* > Am;				// coarse grid matrices for 0 <= m <= amg_levels-1
	std::vector< interpolation_type > Wm;		// interpolation weights for 1 <= m <= amg_levels-1

	// prolongation from m to m-1, m > 0
	void prolongation(int m_from, const std::vector<double>& e, std::vector<double>& v_out) const;
	// restriction from m to m+1
	void restriction(int m_from, const std::vector<double>& e, std::vector<double>& v_out) const;
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
	AMG_Preconditioner(const CSRMatrix* _pA, const SolverParameters& params);
	~AMG_Preconditioner();

	bool SetupPreconditioner() override;
	bool PreconditionedSolve(const std::vector<double>& rhs, std::vector<double>& x) override;
};

#else

class AMG_Preconditioner : public Preconditioner
{
	Preconditioner* Next = nullptr;
public:
	AMG_Preconditioner(const CSRMatrix* _pA, const SolverParameters& params);
	AMG_Preconditioner(const AMG_Preconditioner* prev, int amg_levels, int depth);
	~AMG_Preconditioner();

	bool SetupPreconditioner() override;
	bool PreconditionedSolve(const std::vector<double>& rhs, std::vector<double>& x) override;
};

#endif

class JACOBI_Preconditioner : public Preconditioner
{
public:
	JACOBI_Preconditioner(const CSRMatrix* _pA, const SolverParameters& params) : Preconditioner(_pA, params) {}

	bool SetupPreconditioner() {return true;}
	bool PreconditionedSolve(const std::vector<double>& rhs, std::vector<double>& x) override;
};

class GS_Preconditioner : public Preconditioner
{
public:
	GS_Preconditioner(const CSRMatrix* _pA, const SolverParameters& params) : Preconditioner(_pA, params) {}

	bool SetupPreconditioner() {return true;}
	bool PreconditionedSolve(const std::vector<double>& rhs, std::vector<double>& x) override;
};

class SSOR_Preconditioner : public Preconditioner
{
public:
	SSOR_Preconditioner(const CSRMatrix* _pA, const SolverParameters& params) : Preconditioner(_pA, params) {}

	bool SetupPreconditioner() {return true;}
	bool PreconditionedSolve(const std::vector<double>& rhs, std::vector<double>& x) override;
};

////////////////////////////////////////////////////////////////////

static Preconditioner * CreatePreconditioner(const CSRMatrix* pA, const SolverParameters& parameters)
{
	PreconditionerType ptype = GetPreconditionerType(parameters.GetStringParameter("preconditioner").second);
	if(ptype == PreconditionerType::NONE)
		return new NOOP_Preconditioner(pA, parameters);
	else if(ptype == PreconditionerType::ILU0)
		return new ILU0_Preconditioner(pA, parameters);
	else if(ptype == PreconditionerType::ILUC)
		return new ILUC_Preconditioner(pA, parameters);
	else if(ptype == PreconditionerType::AMG)
		return new AMG_Preconditioner(pA, parameters);
	else if(ptype == PreconditionerType::JACOBI)
		return new JACOBI_Preconditioner(pA, parameters);
	else if(ptype == PreconditionerType::GAUSS_SEIDEL)
		return new GS_Preconditioner(pA, parameters);
	else if(ptype == PreconditionerType::SSOR)
		return new SSOR_Preconditioner(pA, parameters);
	else
	{
		std::cerr << "Sorry, preconditioner was not implemented " << std::endl;
		return NULL;
	}
}


#endif
