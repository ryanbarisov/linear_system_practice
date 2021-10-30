#include <matrix.h>
#include <matrixreader.h>
#include <method.h>
#include <cstring>


double Timer();

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

PreconditionerType GetPreconditionerType(std::string name)
{
	if(name == "AMG")
		return PreconditionerType::AMG;
	else if(name == "ILUC")
		return PreconditionerType::ILUC;
	else
	{
		std::cerr << "Failed to get preconditioner type from string: " << name << std::endl;
		return PreconditionerType::NONE;
	}
}

int main(int argc, char ** argv)
{
	if(argc < 2)
	{
		std::cout << "Usage: ./main <matrix.mtx> <preconditioner_type=AMG|ILUC> [rhs.mtx]" << std::endl;
	}
	else
	{
		std::string filename(argv[1]);
		double t_read = Timer();
		SparseMatrix * A = ReadMatrix(GetMatrixFormat(filename), filename.c_str());
		std::cout << "Time to read matrix " << argv[1] << " " << (Timer() - t_read) << std::endl;
		int n = A->Size();
		std::vector<double> x,b;
		x.resize(n);
		b.resize(n);
		for(int i = 0; i < n; i++)
		{
			b[i] = i;
			x[i] = 0.0;
		}
		if(argc > 3)
		{
			std::cout << "Try to read RHS from MTX file: " << argv[3] << std::endl;
			std::cout << "Finished " << (read_rhs_from_mtx(b, argv[3]) ? "successfully " : "with error.") << std::endl;
		}

		PreconditionerType ptype;
		if(argc > 2)
		{
			std::cout << argv[2] << " is AMG? " << !strcmp(argv[2], "AMG") << " is ILUC? " << !strcmp(argv[2], "ILUC") << std::endl;
			ptype = GetPreconditionerType(std::string(argv[2]));
		}
		bool preconditioned = ptype != PreconditionerType::NONE;
		if(!preconditioned)	std::cout << "Use CG without preconditioning" << std::endl;

		Method * method;
		if(preconditioned) 
			method = new PCG_method();
		else
		 	method = new CG_method();
		method->SetIntegerParameter("maximum_iterations", 100);
		method->SetIntegerParameter("amg_levels", 3);
		method->SetRealParameter("drop_tolerance", 1e-4);
		method->SetIntegerParameter("level_of_fill", 80);

		double t_precond = Timer(), t_solve;
		bool precond = method->Setup(A, ptype);
		t_precond = Timer() - t_precond;

		if(precond)
		{
			std::cout << "Preconditioner time: " << t_precond << std::endl;
			t_solve = Timer();
			bool solve = method->Solve(b,x);
			t_solve = Timer() - t_solve;
			if(solve)
			{
				std::cout << "Solve time: " << t_solve << std::endl;
				SaveVector(x, "solution.txt");
			}
		}

		delete method;
		delete A;
	}

	
	return 0;
}