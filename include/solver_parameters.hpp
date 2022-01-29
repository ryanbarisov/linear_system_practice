#ifndef SOLVER_PARAMETERS_H
#define SOLVER_PARAMETERS_H

#include <map>
#include <string>
#include <utility>


class SolverParameters
{
private:
	static const bool debug = false;
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

#endif // SOLVER_PARAMETERS_H