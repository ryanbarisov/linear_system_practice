#include <matrix.h>
#include <matrixreader.h>
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




int main(int argc, char ** argv)
{
	if(argc < 2)
	{
		//std::cout << "Usage: ./main <matrix.mtx> <preconditioner_type=AMG|ILUC> [rhs.mtx]" << std::endl;
		std::cout << "Usage: ./main <options.json>" << std::endl;
	}
	else if(true)
	{
		std::string filename(argv[1]);

		rapidjson::Document options = parseJsonFromFile(filename);
		if(options.HasParseError())
		{
			std::cerr << "Error while opening options file: " << filename 
				<< ", error code (see rapidjson/error/error.h): " << options.GetParseError() << std::endl;
			return 1;
		}
		assert(options.HasMember("matrix"));

		std::string matrix_filename = options["matrix"].GetString();
		double t_read = Timer();
		SparseMatrix * A = ReadMatrix(GetMatrixFormat(matrix_filename), matrix_filename.c_str());
		std::cout << "Time to read matrix " << matrix_filename << " " << (Timer() - t_read) << std::endl;

		int n = A->Size();
		std::vector<double> x,b;
		x.resize(n);
		b.resize(n);
		std::fill(x.begin(), x.end(), 0.0);
		if(options.HasMember("rhs"))
		{
			std::string rhs_filename = options["rhs"].GetString();
			std::cout << "Try to read RHS from MTX file: " << rhs_filename << std::endl;
			std::cout << "Finished " << (read_rhs_from_mtx(b, rhs_filename.c_str()) ? "successfully " : "with error.") << std::endl;
		}
		else
		{
			for(int i = 0; i < n; i++) b[i] = i;
		}

		PreconditionerType ptype = PreconditionerType::NONE;
		if(options.HasMember("preconditioner"))
		{
			ptype = GetPreconditionerType(options["preconditioner"].GetString());
		}
		bool preconditioned = ptype != PreconditionerType::NONE;
		if(!preconditioned)	std::cout << "Use CG without preconditioning" << std::endl;

		Method * method;
		if(preconditioned) 
			method = new PCG_method();
		else
			method = new CG_method();

		if(options.HasMember("parameters"))
		{
			const rapidjson::Value& parameters = options["parameters"];
			if(parameters.HasMember("maximum_iterations"))
				method->SetIntegerParameter("maximum_iterations", parameters["maximum_iterations"].GetInt());
			// ILUC
			if(parameters.HasMember("level_of_fill"))
				method->SetIntegerParameter("level_of_fill", parameters["level_of_fill"].GetInt());
			if(parameters.HasMember("drop_tolerance"))
				method->SetRealParameter("drop_tolerance", parameters["drop_tolerance"].GetDouble());
			// AMG
			if(parameters.HasMember("amg_levels"))
				method->SetIntegerParameter("amg_levels", parameters["amg_levels"].GetInt());
			if(parameters.HasMember("niters_smooth"))
				method->SetIntegerParameter("niters_smooth", parameters["niters_smooth"].GetInt());
			if(parameters.HasMember("w_cycle"))
				method->SetIntegerParameter("w_cycle", parameters["w_cycle"].GetInt());
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
			}
		}

		delete method;
		delete A;
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
		method->SetIntegerParameter("niters_smooth", 5);
		method->SetIntegerParameter("w_cycle", 0);
		method->SetRealParameter("drop_tolerance", 1.0e-4);
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