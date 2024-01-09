#ifndef MATRIX_READER_H
#define MATRIX_READER_H

#include <matrix.h>


class MatrixReader
{
public:
	static CSRMatrix* ReadMatrix(const char * filename);
};

class MTXMatrixReader : public MatrixReader
{
public:
	static CSRMatrix* ReadMatrix(const char * filename);
};

class CSRMatrixReader : public MatrixReader
{
public:
	static CSRMatrix* ReadMatrix(const char * filename);
};


CSRMatrix* ReadMatrix(MatrixFormat fmt, const char * filename);

bool read_vector_mtx(std::vector<double>& rhs, const char* filename);


#endif // MATRIX_READER_H
