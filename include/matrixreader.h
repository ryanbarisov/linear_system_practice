#ifndef MATRIX_READER_H
#define MATRIX_READER_H

#include <matrix.h>


class MatrixReader
{
public:
	static SparseMatrix* ReadMatrix(const char * filename);
};

class MTXMatrixReader : public MatrixReader
{
public:
	static SparseMatrix* ReadMatrix(const char * filename);
};

class CSRMatrixReader : public MatrixReader
{
public:
	static SparseMatrix* ReadMatrix(const char * filename);
};


SparseMatrix* ReadMatrix(MatrixFormat fmt, const char * filename);

bool read_rhs_from_mtx(std::vector<double>& rhs, const char* filename);


#endif // MATRIX_READER_H
