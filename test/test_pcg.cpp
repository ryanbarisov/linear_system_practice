#include <matrix.h>
#include <matrixreader.h>
#include <matrixwriter.h>
#include <method.h>
#include <cstring>

double Timer();


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
	bool print_info = argc == 1;
	std::string matrix_filename, rhs_filename;
	SolverParameters params;
	for(int i = 1; i < argc && !print_info; i+=2)
	{
		if(strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "--matrix") == 0)
			matrix_filename = argv[i+1];
		else if(strcmp(argv[i], "-b") == 0 || strcmp(argv[i], "--rhs") == 0)
			rhs_filename = argv[i+1];
		else if(strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--params") == 0)
			params.Setup(argv[i+1]);
		else print_info = true;
	}
	if(print_info)
	{
		std::cout << "Usage: " << argv[0] << " -m matrix.mtx -p params.txt [-b rhs.mtx]" << std::endl;
		return 0;
	}

	SparseMatrix * A = ReadMatrix(GetMatrixFormat(matrix_filename), matrix_filename.c_str());

	int n = A->Size();
	std::vector<double> x(n,0.0),b(n);
	read_vector_mtx(b, rhs_filename.c_str());

	Method * method = CreateMethod(params);
	if(!method) return 1;
	double t_prec = Timer();
	bool precond = method->Setup(A, params);
	t_prec = Timer()-t_prec;

	if(precond)
	{
		std::cout << "Preconditioner time: " << t_prec << std::endl;
		double t_solve = Timer();
		bool solve = method->Solve(b,x);
		t_solve = Timer() - t_solve;
		if(solve)
		{
			std::cout << "Solve time: " << t_solve << std::endl;
			SaveVector(x, "solution.txt");
			//if(file_solution) print_norms(x, x_file);
		}
	}

	delete method;
	delete A;
	return 0;
}
