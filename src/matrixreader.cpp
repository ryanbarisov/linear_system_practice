#include <matrix.h>
#include <matrixreader.h>
#include <matrixwriter.h>

SparseMatrix* ReadMatrix(MatrixFormat fmt, const char * filename)
{
	if(fmt == MatrixFormat::MTX)
		return MTXMatrixReader::ReadMatrix(filename);
	else if(fmt == MatrixFormat::CSR)
		return CSRMatrixReader::ReadMatrix(filename);
	else
	{
		std::cerr << "Matrix format is not supported." << std::endl;
		return nullptr;
	}
}

SparseMatrix* MTXMatrixReader::ReadMatrix(const char * filename)
{
	SparseMatrix* A = nullptr;

	std::ifstream ifs(filename);
	if(!ifs.is_open())
	{
		std::cout << "Failed to open file " << filename << std::endl;
	}
	else
	{
		int rows, cols, nnz,  i, j, k = 0;
		bool read_first_line = false;
		double aij;
		std::string s;
		while(!ifs.eof())
		{
			std::getline(ifs,s);
			if(s.empty())	continue;
			else if(s[0] == '%')
				continue;	// skip comment line
			else
			{
				std::stringstream ss(s);
				if(!read_first_line)
				{
					ss >> rows >> cols >> nnz;
					A = new SparseMatrix(rows);
				}
				else
				{
					ss >> i >> j >> aij;
					//(*A)[i-1].add_element(j-1, aij); k++;
					{
						// for symmetric matrices from Florida collection
						(*A)[i-1].add_element(j-1, aij); k++;
						(*A)[j-1].add_element(i-1, aij); k++;
					}
				}
				if(ss.fail())
				{
					std::cout << "Error occurred during reading string: " << s << std::endl;
					break;
				}
				if(!read_first_line)
					read_first_line = true;
			}
		}
		std::cout << "Matrix " << filename << ": " << rows << " X " << cols << ", nnz=" << k << std::endl;
	}
	ifs.close();

	return A;
}




SparseMatrix* CSRMatrixReader::ReadMatrix(const char * filename)
{
	SparseMatrix* A = nullptr;

	std::ifstream ifs(filename);
	if(!ifs.is_open())
	{
		std::cout << "Failed to open file " << filename << std::endl;
	}
	else
	{
		int rows, nnz;
		ifs >> rows;

		int * ia = new int[rows+1];
		for(int i = 0; i < rows+1; i++)
			ifs >> ia[i];

		nnz = ia[rows] - ia[0];
		int * ja = new int[nnz];
		double * a = new double[nnz];
		for(int i = 0; i < nnz; i++) ifs >> ja[i];
		for(int i = 0; i < nnz; i++) ifs >> a[i];
		
		A = new SparseMatrix(rows);
		for(int r = 0; r < rows; r++)
		{
			sparse_row& row = (*A)[r];
			for(int k = ia[r]-1; k < ia[r+1]-1; k++)
				row.add_element(ja[k]-1, a[k]);
		}

		delete[] ia;
		delete[] ja;
		delete[] a;
	}
	ifs.close();

	return A;
}




bool read_rhs_from_mtx(std::vector<double>& rhs, const char* filename)
{
	bool success = false;
	std::ifstream ifs(filename);
	if(!ifs.is_open())
	{
		std::cout << "Failed to open file " << filename << std::endl;
	}
	else
	{
		bool flag = false;
		int i, j, k = 0, M,N,L;
		double aij;
		std::string s;
		while(!ifs.eof())
		{
			std::getline(ifs,s);
			if(s.empty())	continue;
			else if(s[0] == '%') continue;	// skip comment 
			else
			{
				std::stringstream ss(s);
				if(!flag)
				{
					ss >> M >> N;
					assert(N == 1);
				}
				else
				{
					ss >> aij;
					rhs[k++] = aij;
				}
				if(ss.fail())
				{
					std::cout << "Error occurred during reading string: " << s << std::endl;
					break;
				}
				if(!flag)
				{
					rhs.resize(M);
					flag = true;
				}
			}
		}
		success = k == M;
	}
	ifs.close();

	return success;
}