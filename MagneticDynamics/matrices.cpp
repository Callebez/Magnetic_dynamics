#include "matrices.hpp"

matrix::matrix(uint nrows, uint ncolumns)
{
    double ** m = new double*[rows];
    if(rows)
    {
        m[0] = new double[nrows*ncolumns];
        for (uint i = 1; i < nrows; ++i)
        {
            m[i] = m[0] + i * ncolumns;
        }
    }
    data = m; 
    rows = nrows;
    columns = nrows;
}
unsigned int matrix::get_rows()
{
	return rows;
}
unsigned int matrix::get_columns()
{
	return columns;
}
double matrix::get_value(unsigned int row_position, unsigned int column_position)
{
	if (row_position > rows || column_position > columns)
	{
		throw std::invalid_argument("Invalid row/column entrance!");
	}
	return data[row_position][column_position];
}
void matrix::set_rows(unsigned int new_nrows)
{
	if (new_nrows == 0)
	{
		throw std::invalid_argument("number of rows cannot be zero!");
	}
	rows = new_nrows;
}
void matrix::set_columns(unsigned int new_ncolumns)
{
	if (new_ncolumns == 0)
	{
		throw std::invalid_argument("number of columns cannot be zero!");
	}
	columns = new_ncolumns;
}
void matrix::print_matrix()
{
	if (rows == 0 || columns == 0)
	{
		throw std::invalid_argument("Invalid matrix!");
	}
	for (unsigned int i = 0; i < get_rows(); i++)
	{
		for (unsigned int j = 0; j < get_columns(); j++)
		{
			std::cout << matrix::get_value(i,j) << " ";
		}
		std::cout << std::endl;
	}
}
std::vector<double> matrix::apply_matrix(matrix& A, std::vector<double>& v)
{
	if (A.get_columns() != v.size())
	{
		throw std::invalid_argument("The matrix A and the vector v must have the same dimension!");
	}
	double sum = 0;

	std::vector<double> res (v.size());
	for (unsigned int i = 0; i < A.get_rows(); i++)
	{
		sum = 0;
		for (unsigned int j = 0; j < A.get_columns(); j++)
		{
			sum += A.get_value(i,j) * v[j];
		}
		res[i] = sum;
	}
	return res;
}

void matrix::set_value(unsigned int row_position, unsigned int column_position, double value)
{
	if (row_position > rows)
	{
		throw std::invalid_argument("Invalid row position");
	}
	if (column_position > columns)
	{
		throw std::invalid_argument("Invalid column position");
	}
	data[row_position][column_position] = value;
}
matrix matrix::sum_matrix(matrix& A, matrix& B)
{
	if (A.get_columns() != B.get_columns())
	{
		throw std::invalid_argument("matrices must have the same number of columns");
	}
	if (A.get_rows() != B.get_rows())
	{
		throw std::invalid_argument("matrices must have the same number of rows");
	}
	matrix C(A.get_rows(), A.get_columns());
	for (unsigned int i = 0; i < A.get_rows(); i++)
	{
		for (unsigned int j = 0; j < A.get_columns(); j++)
		{
			C.set_value(i, j, A.get_value(i, j) + B.get_value(i, j));
		}
	}
	return C;
}
matrix matrix::sub_matrix(matrix& A, matrix& B)
{
	if (A.get_columns() != B.get_columns())
	{
		throw std::invalid_argument("matrices must have the same number of columns");
	}
	if (A.get_rows() != B.get_rows())
	{
		throw std::invalid_argument("matrices must have the same number of rows");
	}
	matrix C(A.get_rows(), A.get_columns());
	for (unsigned int i = 0; i < A.get_rows(); i++)
	{
		for (unsigned int j = 0; j < A.get_columns(); j++)
		{
			C.set_value(i, j, A.get_value(i, j) - B.get_value(i, j));
		}
	}
	return C;
}
matrix matrix::mult_matrix(matrix& A, matrix& B)
{
	matrix C(A.get_rows(), B.get_columns());
	double sum = 0;
	if (A.get_columns() != B.get_rows())
	{	
		std::cout << A.get_columns();
		std::cout << B.get_rows();
		throw std::invalid_argument("number of columns of A must be the same as the number of rows in B.");
	}
	for (unsigned int i = 0; i < C.get_rows(); ++i)
	{
		for (unsigned int j = 0; j < C.get_columns(); ++j)
		{
			sum = 0;
			for (unsigned int k = 0; k < A.get_columns(); ++k)
			{
				sum += A.get_value(i,k) * B.get_value(k,j);
			}
			C.set_value(i,j,sum);
		}
	}
	return C;
}
void matrix::print_vector(std::vector<double>& v)
{
	std::cout << "(" << v[0];
	for (unsigned int i = 1; i < v.size(); i++)
	{
		std::cout << ", " << i;
	}
	std::cout << ")\n";
}
std::vector<double> matrix::sum_vectors(std::vector<double>& v, std::vector<double>u)
{
	if (v.size() != u.size())
	{
		throw std::invalid_argument("Vectors must have the same dimension!");
	}
	std::vector<double> res(v.size());
	for (unsigned int i = 0; i < v.size(); i++)
	{
		res[i] = v[i] + u[i];
	}
	return res;
}
std::vector<double> matrix::sub_vectors(std::vector<double>& v, std::vector<double>u)
{
	if (v.size() != u.size())
	{
		throw std::invalid_argument("Vectors must have the same dimension!");
	}
	std::vector<double> res(v.size());
	for (unsigned int i = 0; i < v.size(); i++)
	{
		res[i] = v[i] - u[i];
	}
	return res;
}
// Euclidian inner product
double matrix::dotProduct(std::vector<double>& v, std::vector<double>& u)
{
	if (v.size() != u.size())
	{
		throw std::invalid_argument("Vectors must have the same dimension!");
	}
	double prod = 0;
	for (unsigned int i = 0; i < v.size(); i++)
	{
		prod += v[i] * u[i];
	}
	return prod;
}
double matrix::norm(std::vector<double>& v){ return sqrt(dotProduct(v, v)); }
double matrix::distance(std::vector<double>& v, std::vector<double>& u) 
{
	std::vector<double> a = sub_vectors(v, u);
	return norm(a);
}
// double* matrix::axpy(double*& x, double*& y, double a, uint size )
// {
// 	double* res = new double[size];
// 	for (unsigned int i = 0; i < size; i++)
// 	{
// 		res[i] = a * x[i] + y[i];
// 	}
// 	return res;
// }
void matrix::axpy(double*& x, double*& y, double a, uint size,double*& res )
{	
	for (unsigned int i = 0; i < size; i++)
	{
		res[i] = a * x[i] + y[i];
	}
}
std::vector<double> matrix::ax(std::vector<double>& x, double a)
{
	std::vector<double> res(x.size());
	for (unsigned int i = 0; i < x.size(); i++)
	{
		res[i] = a * x[i];
	}
	return res;
}

matrix matrix::copy_matrix()
{
	matrix v = matrix(rows, columns);
	for (unsigned int i = 0; i < rows; i++)
	{
		for (unsigned int j = 0; j < columns; j++)
		{
			v.data[i][j] = data[i][j];
		}
	}
	return v;
}

matrix matrix::identityMatrix(unsigned int order)
{
	matrix id = matrix(order, order);
	for (unsigned int i = 0; i < order; i++)
	{
		for (unsigned int j = 0; j < order; j++)
		{
			id.data[i][j] = 0;
		}
		id.data[i][i] = 1;
	}
	return id;
}
void matrix::transpose()
{
	matrix aux = copy_matrix();
	for (unsigned int i = 0; i < rows; i++)
	{
		for (unsigned int j = 0; j < columns; j++)
		{
			data[i][j] = aux.data[j][i];
		}
	}
}
void matrix :: printMatrixToFile(std::vector<std::pair<std::vector<double>, double>>& matrix, std::string fileName)
{
	std::ofstream output;

	output.open(fileName);
	for (auto i : matrix)
	{
		output << i.second;
		for (unsigned int j = 0; j < i.first.size(); j++)
		{
			output << " " << i.first[j];
		}
		output << std::endl;
	}
	output.close();
}