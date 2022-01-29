#include <matrix.h>
#include <matrixreader.h>
#include <matrixwriter.h>
#include <method.h>
#include <cstring>
#include <rapidjson/rapidjson.h>
#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>

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
	else if(name == "ILU0")
		return PreconditionerType::ILU0;
	else if(name == "JAC")
		return PreconditionerType::JACOBI;
	else if(name == "GS")
		return PreconditionerType::GAUSS_SEIDEL;
	else
	{
		std::cerr << "Failed to get preconditioner type from string: " << name << std::endl;
		return PreconditionerType::NONE;
	}
}



static rapidjson::Document parseJsonFromFile(std::string filename)
{
	rapidjson::Document options;
	std::ifstream ifs(filename.c_str());
	if(!ifs.good())
	{
		std::cerr << "Failed to open options file: " << filename << std::endl;
	}
	else
	{
		rapidjson::IStreamWrapper isw(ifs);
		options.ParseStream(isw);
		ifs.close();
	}
	return options;
}


SparseMatrix * HilbertMatrix(int n)
{
	SparseMatrix * A = new SparseMatrix(n);
	for(int i = 0; i < n; i++)
	{
		sparse_row& r = (*A)[i];
		for(int j = 0; j < n; j++)
		{
			r.add_element(j, 1.0/(i+j+1+1));
		}
	}

	return A;
}

SparseMatrix * LaplaceMatrix(int n)
{
	SparseMatrix * A = new SparseMatrix(n);
	for(int i = 0; i < n; i++)
	{
		sparse_row& r = (*A)[i];
		for(int j = 0; j < 3; j++)
		{
			if(i-1 >= 0) r.add_element(i-1, -1.0);
			r.add_element(i, 2.0);
			if(i+1 < n)  r.add_element(i+1, -1.0);
		}
	}

	return A;
}

bool read_vector(const rapidjson::Document& options, std::string str, std::vector<double>& v, int n)
{
	bool read_v = false;
	if(options.HasMember(str.c_str()))
	{
		std::string rhs_filename = options[str.c_str()].GetString();
		std::cout << "Try to read RHS from MTX file: " << rhs_filename << std::endl;
		read_v = read_rhs_from_mtx(v, rhs_filename.c_str()); 
		std::cout << "Finished " << (read_v ? "successfully " : "with error.") << std::endl;
	}
	if(!read_v)
	{
		std::cout << "Provide analytical vector" << std::endl;
		double pi = 3.1415926535;
		for(int i = 0; i < n; i++) v[i] = sin(pi / n * i)/n * (1.0*rand()/RAND_MAX);
	}
	return read_v;
}

void print_norms(const std::vector<double>& x1, const std::vector<double>& x2)
{
	int n1 = x1.size(), n2 = x2.size();
	assert(n1 == n2);
	int n = std::min(n1, n2);
	double l1 = 0, l2 = 0, linf = 0;
	for(int i = 0; i < n; i++)
	{
		double diff = fabs(x1[i]-x2[i]);
		l1 += diff;
		l2 += diff*diff;
		linf = std::max(diff, linf);
	}
	l1 /= n;
	l2 = sqrt(l2)/n;
	std::cout << "Norms: L2 " << l2 << "; Linf " << linf << "; L1 " << l1 << std::endl;
}



int main(int argc, char ** argv)
{
	if(argc < 2)
	{
		std::cout << "Usage: ./main <options.json>" << std::endl;
		return 0;
	}
	
	std::string filename(argv[1]);

	rapidjson::Document options = parseJsonFromFile(filename);
	if(options.HasParseError())
	{
		std::cerr << "Error while opening options file: " << filename 
			<< ", error code (see rapidjson/error/error.h): " << options.GetParseError() << std::endl;
		return 1;
	}
	SparseMatrix * A;
	if(options.HasMember("matrix"))
	{
		std::string matrix_filename = options["matrix"].GetString();
		double t_read = Timer();
		A = ReadMatrix(GetMatrixFormat(matrix_filename), matrix_filename.c_str());
		std::cout << "Time to read matrix " << matrix_filename << " " << (Timer() - t_read) << std::endl;
	}
	else
	{
		int n = 40;
		std::cerr << "No matrix provided, use Laplace matrix of size " << n << std::endl;
		A = LaplaceMatrix(n);
		//MTXMatrixWriter::WriteMatrix(*A, "laplace.mtx");
	}
	

	int n = A->Size();
	std::vector<double> x(n,0.0),b(n);
	std::vector<double> x_file(n);
	read_vector(options, "rhs", b, n);
	bool file_solution = read_vector(options, "sol", x_file, n);

	PreconditionerType ptype = PreconditionerType::NONE;
	if(options.HasMember("preconditioner"))
	{
		ptype = GetPreconditionerType(options["preconditioner"].GetString());
		std::cout << "Use preconditioner " << options["preconditioner"].GetString() << std::endl;
	}
	bool preconditioned = ptype != PreconditionerType::NONE;
	if(!preconditioned)	std::cout << "Use CG without preconditioning" << std::endl;

	Method * method;
	if(options.HasMember("solver"))
	{
		std::string solver_type = options["solver"].GetString();
		if(solver_type == "CG")
		{
			if(preconditioned) method = new PCG_method();
			else method = new CG_method();
		}
		else
		{
			if(preconditioned) method = new PBICGStab_method();
			else method = new BICGStab_method();
		}
	}
	

	if(options.HasMember("parameters"))
	{
		const rapidjson::Value& parameters = options["parameters"];
		if(parameters.HasMember("maximum_iterations"))
			method->SetIntegerParameter("maximum_iterations", parameters["maximum_iterations"].GetInt());
		if(ptype == PreconditionerType::ILUC)
		{
			if(parameters.HasMember("level_of_fill"))
				method->SetIntegerParameter("level_of_fill", parameters["level_of_fill"].GetInt());
			if(parameters.HasMember("drop_tolerance"))
				method->SetRealParameter("drop_tolerance", parameters["drop_tolerance"].GetDouble());
		}
		if(ptype == PreconditionerType::AMG)
		{
			if(parameters.HasMember("amg_levels"))
				method->SetIntegerParameter("amg_levels", parameters["amg_levels"].GetInt());
			if(parameters.HasMember("niters_smooth"))
				method->SetIntegerParameter("niters_smooth", parameters["niters_smooth"].GetInt());
			if(parameters.HasMember("w_cycle"))
				method->SetIntegerParameter("w_cycle", parameters["w_cycle"].GetInt());
		}
	}

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
			if(file_solution) print_norms(x, x_file);
		}
	}

	delete method;
	delete A;
}
