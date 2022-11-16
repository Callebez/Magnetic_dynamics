#include <iostream>
#include<fstream>
#include "matrices.h"
#include "functions.h"
void plot3D(const char* fileName, const char* inputFile, const char* systemName, const char* keys)
{
	FILE* gnupipe = _popen("gnuplot -persist", "w");
	if (gnupipe)
	{
		fprintf_s(gnupipe,"set terminal pngcairo enhanced size 1080,720 font \'Heveltica, 15\ '\n\n");
		fprintf_s(gnupipe, "set title \"%s \" font \'Helvetica, 16\' \n\n ", systemName);
		fprintf_s(gnupipe, "set tics font \'Helvetica, 14\' \n\n");
		fprintf_s(gnupipe, "set border lw 3 \n \n");
		fprintf_s(gnupipe, "set key opaque \n\n");
		fprintf_s(gnupipe,"set output \'%s.png\'\n", fileName);
		fprintf_s(gnupipe, "splot \'%s.txt\' u 2:3:4 w lines lw 0.5 lc 7 title \'%s\'\n", fileName, keys);
	}

}
void printMatrixToFile(std::vector<std::pair<std::vector<double>, double>>& matrix, std::string fileName)
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
	std::cout << fileName;
}

int main(void)
{

	std::vector<double> params = {9.8,1.0, 100.0};

	function f(function::sistemaAmigos, params, 2, 2);
	std::vector<double> x = { 10.0,0.000};
	
	std::vector<std::pair<std::vector<double>, double>> y = function::runge_kutta4th(f, make_pair(x, 0), 100, 0.01);

	printMatrixToFile(y, "amigos.txt");

	return 0;
}