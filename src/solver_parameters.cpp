#include <solver_parameters.h>
#include <fstream>
#include <sstream>
#include <cstdlib>


// check: 0 -- integer, 1 -- double
static bool is_number(const std::string& s, int check = 0)
{
	if(s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+'))) return false;

	char * p;
	if(check == 0)
		strtol(s.c_str(), &p, 10);
	else
		strtod(s.c_str(), &p);

	return (*p == 0);
}

void SolverParameters::SetParameter(std::string name, std::string value)
{
	if(is_number(value, 0)) // integer?
		SetIntegerParameter(name, std::atoi(value.c_str()));
	else if(is_number(value, 1)) // double?
		SetRealParameter(name, std::atof(value.c_str()));
	else
		SetStringParameter(name, value);
}

// parameters are given as key-value pairs separated by whitespaces
// lines starting with # are comment lines
void SolverParameters::Setup(std::string filename)
{
	SetupDefault();
	// default parameters
	std::ifstream ifs(filename);
	while(ifs.good())
	{
		std::string line;
		std::getline(ifs, line);
		if(line[0] == '#') continue;
		std::stringstream ss(line);
		std::string key, value;
		ss >> key >> value;
		SetParameter(key, value);
	}
}

void SolverParameters::SetupDefault()
{
	SetStringParameter("preconditioner", "none");
	SetStringParameter("solver", "cg");
	// solver parameters
	SetIntegerParameter("maximum_iterations", 1000);
	SetRealParameter("relative_tolerance", 1.0e-8);
	SetRealParameter("absolute_tolerance", 1.0e-12);
	// default preconditioner parameters
	// ILUC
	SetRealParameter("drop_tolerance", 1.0e-4);
	SetIntegerParameter("level_of_fill", 80);
	// AMG
	SetIntegerParameter("amg_levels", 2);
	SetIntegerParameter("niters_smooth", 2);
	SetIntegerParameter("w_cycle", 0);
}
