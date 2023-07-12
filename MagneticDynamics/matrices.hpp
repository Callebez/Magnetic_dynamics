#pragma once
#include <vector>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <complex>
#include <cmath>
#include <memory>
#include <algorithm>
#include <functional>
class matrix
{
	private:
		unsigned int rows = 0, columns = 0;
		// std::vector<std::vector<double>> data;

	public:
		double** data;

		matrix(unsigned int nrows, unsigned int ncolumns);
		
		unsigned int get_rows();
		unsigned int get_columns();
		double** get_data(){return data;};
		double get_value(unsigned int row_position, unsigned int column_position);
		
		void set_rows(unsigned int new_nrows);
		void set_columns(unsigned int new_columns);
		void set_value(unsigned int row_position, unsigned int column_position, double value);
		
		void print_matrix();
		matrix copy_matrix();
		void transpose();
		static matrix sum_matrix(matrix& A, matrix& B);
		static matrix sub_matrix(matrix& A, matrix& B);
		static matrix mult_matrix(matrix& A, matrix& B);
		static std::vector<double> apply_matrix(matrix& A, std::vector<double>& v);

		static void print_vector(std::vector<double>& v);
		static double distance(std::vector<double>& v, std::vector<double>& u);
		static double norm(std::vector<double>& v);
		static double dotProduct(std::vector<double>& v, std::vector<double>& u);
		static std::vector<double> sub_vectors(std::vector<double>& v, std::vector<double>u);
		static std::vector<double> sum_vectors(std::vector<double>& v, std::vector<double>u);
		static void axpy(double*& x,double*& y, double a, uint size,double*& res);
		static std::vector<double> ax(std::vector<double>& x, double a);
		static matrix identityMatrix(unsigned int order);
		static void printMatrixToFile(std::vector<std::pair<std::vector<double>, double>>& matrix, std::string fileName);

	};

