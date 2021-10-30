#ifndef MATRIX_WRITER_H
#define MATRIX_WRITER_H

#include <matrix.h>


class MatrixWriter
{
public:
	static bool WriteMatrix(const SparseMatrix& mat, const char * filename);
};

class MTXMatrixWriter : public MatrixWriter
{
public:
	static bool WriteMatrix(const SparseMatrix& mat, const char * filename);
};

class CSRMatrixWriter : public MatrixWriter
{
public:
	static bool WriteMatrix(const SparseMatrix& mat, const char * filename);
};


bool WriteMatrix(const SparseMatrix& mat, const char * filename, MatrixFormat fmt);


#endif // MATRIX_WRITER_H