#include <matrix.h>
#include <matrixwriter.h>


bool WriteMatrix(const SparseMatrix& mat, const char * filename, MatrixFormat fmt)
{
	if(fmt == MatrixFormat::MTX)
		return MTXMatrixWriter::WriteMatrix(mat, filename);
	//else if(fmt == MatrixFormat::CSR)
	//	return CSRMatrixWriter::WriteMatrix(mat, filename);
	else
	{
		std::cerr << "Matrix format is not supported." << std::endl;
		return false;
	}
}

bool MTXMatrixWriter::WriteMatrix(const SparseMatrix& mat, const char * filename)
{
	std::ofstream ofs(filename);
	ofs << "%% MTX Format matrix" << std::endl;
	ofs << "%% Linear System Practice!" << std::endl;
	for(int i = 0; i < mat.Size(); i++)
	{
		const sparse_row& r = mat[i];
		for(int k = 0; k < r.row.size(); k++)
		{
			int col = r.row[k].first;
			double val = r.row[k].second;
			ofs << i+1 << " " << col+1 << " " << val << std::endl;
		}
	}
	ofs.close();
	return true;
}