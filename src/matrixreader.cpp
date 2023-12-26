#include <matrix.h>
#include <matrixreader.h>
#include <matrixwriter.h>
#include <algorithm>
#include <cctype>

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

static bool my_isspace(char ch)
{
    return std::isspace(static_cast<unsigned char>(ch));
}

// split string on words using space as separator
// return amount of words
static int parse_string_numbers(std::string str)
{
	// doesn't work correctly if words separated by several spaces
	//return std::count(str.begin(), str.end(), ' '); 
	// naive approach
	std::string::iterator beg = str.begin(), end = str.end();
	std::string::iterator it;
	int count = 0;
	bool word = false, space = false;
	// every non-space character is a word character 
	for(it = beg; ; ++it)
	{
		if(it == end)
		{
			if(word) count++;
			break;
		}
		space = my_isspace(*it);
		// space finished on previous symbol
		if(!space && !word)
		{
			word = true;
		}
		// word finished on previous symbol, increment word counter
		else if(space && word)
		{
			word = false;
			count++;
		}

		// test redundancy
		if(space && word)
			std::cout << *it << std::endl;
	}
	return count;
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
					std::string line = ss.str();
					int nums = parse_string_numbers(line);
					if(nums == 3)
						ss >> rows >> cols >> nnz;
					else if(nums == 2)
					{
						ss >> rows >> cols;
						nnz = rows;
					}
					else if(nums == 1)
					{
						ss >> rows;
						cols = 1;
						nnz = rows;
					}

					A = new SparseMatrix(rows);
				}
				else
				{
					ss >> i >> j >> aij;
					// (*A)[i-1].add_element(j-1, aij); k++;
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
		//MTXMatrixWriter::WriteMatrix(*A, "save.mtx");
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




bool read_vector_mtx(std::vector<double>& rhs, const char* filename)
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
					int nums = parse_string_numbers(s);
					if(nums == 2)
						ss >> M >> N;
					else if(nums == 1)
					{
						ss >> M;
						N = 1;
					}
					
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
