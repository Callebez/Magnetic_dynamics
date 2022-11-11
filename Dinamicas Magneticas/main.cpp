#include <iostream>
#include "matrices.h"

int main(void)
{
	unsigned int nrows= 2, ncolumns=2;
	matrix M(nrows, ncolumns);
	M.set_value(0, 0, 2);
	M.set_value(0, 1, 1);
	//M.set_value(0, 2, 4);
	M.set_value(1, 0, 0);
	M.set_value(1, 1, 1);
	//M.set_value(1, 2, 1);
	M.print_matrix();
	matrix A = matrix::mult_matrix(M, M);
	A.print_matrix();
	std::vector<double> v = { 0, 1 };
	std::vector<double> a = matrix::apply_matrix(A, v);
	
	matrix::print_vector(a);



	return 0;
}