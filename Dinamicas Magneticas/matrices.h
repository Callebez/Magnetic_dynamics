#pragma once
#include <vector>
#include<cstdint>
#include<iostream>

class matrix
{
	private:
		unsigned int rows, columns;
		double** data;

	public:
		matrix(unsigned int nrows, unsigned int ncolumns);
		
		unsigned int get_rows();
		unsigned int get_columns();
		double get_value(unsigned int row_position, unsigned int column_position);
		
		void set_rows(unsigned int new_nrows);
		void set_columns(unsigned int new_columns);
		void set_value(unsigned int row_position, unsigned int column_position, double value);
		
		void print_matrix();

		static matrix sum_matrix(matrix& A, matrix& B);
		static matrix sub_matrix(matrix& A, matrix& B);
		static matrix mult_matrix(matrix& A, matrix& B);
		static std::vector<double> apply_matrix(matrix& A, std::vector<double>& v);

		static void print_vector(std::vector<double>& v);
		static double distance(std::vector<double>& v, std::vector<double>& u);
		static double norm(std::vector<double>& v);
		static double inner_prod(std::vector<double>& v, std::vector<double>& u);
		static std::vector<double> sub_vectors(std::vector<double>& v, std::vector<double>u);
		static std::vector<double> sum_vectors(std::vector<double>& v, std::vector<double>u);
};

