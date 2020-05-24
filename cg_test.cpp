#include "matrix.hpp"
#include "method.hpp"
#include <cstring>




int main2(int argc, char ** argv)
{
	// if(argc < 2)
	// {
	// 	std::cout << "Usage: ./main <csr_matrix.mat>" << std::endl;
	// }
	// else
	// {
		const int n = 100;
		std::vector<double> a;
		std::vector<int> ia,ja;
		
		ia.resize(n+1);
		int nnz = (n-2)*3 + 2*2;
		a.resize(nnz);
		ja.resize(nnz);

		ia[0] = 1;
		for(int i = 1; i < n+1; i++)
		{
			ia[i] = ia[i-1] + ((i == 1 || i == n) ? 2 : 3);
		}
		int j = 0;
		for(int i = 0; i < n; i++)
		{
			if(i == 0)
			{
				a[j] = 2;
				a[j+1] = -1;
				ja[j] = i+1;
				ja[j+1] = i+2;
				j += 2;
			}
			else if(i == n-1)
			{
				a[j] = -1;
				a[j+1] = 2;
				ja[j] = i;
				ja[j+1] = i+1;
				j += 2;
			}
			else
			{
				a[j] = -1;
				a[j+1] = 2;
				a[j+2] = -1;
				ja[j] = i;
				ja[j+1] = i+1;
				ja[j+2] = i+2;
				j += 3;
			}
		}

		CSR_matrix A(a,ia,ja);
		A.Save("out.dat");
	// }
	return 0;
}

int main3(int argc, char ** argv)
{
	if(argc < 2)
	{
		std::cout << "Usage: ./main <csr_matrix.mat> <preconditioner_type=AMG|ILUC>" << std::endl;
	}
	else
	{
		MTX_matrix * A = new MTX_matrix(argv[1]);
		int n = A->Size();
		std::vector<double> x,b;
		x.resize(n);
		b.resize(n);
		for(int i = 0; i < n; i++)
		{
			b[i] = i;
			x[i] = 0.0;
		}

		PreconditionerType ptype;
		if(argc > 2)
		{
			std::cout << argv[2] << " is AMG? " << !strcmp(argv[2], "AMG") << " is ILUC? " << !strcmp(argv[2], "ILUC") << std::endl;
			if(!strcmp(argv[2], "AMG"))	
				ptype = PreconditionerType::AMG;
			else if(!strcmp(argv[2], "ILUC"))
				ptype = PreconditionerType::ILUC;
			else
				ptype = PreconditionerType::NONE;
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


		if(method->Setup(A, ptype) && method->Solve(b,x))
		{
			SaveVector(x, "solution.txt");
		}
		A->Save("A.mtx");

		delete method;
		delete A;
		std::cin.ignore();
	}

	
	return 0;
}


int main(int argc, char ** argv)
{
	return main3(argc,argv);
}