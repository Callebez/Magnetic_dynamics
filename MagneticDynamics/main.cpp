#include <iostream>
#include<fstream>
#include<chrono>
#include<cstdlib>
#include "matrices.h"
#include "functions.h"
#include <ctime>
#include "plotting.h"

#define PI 3.14159265359

int main(void)
{
	auto start = std::chrono::high_resolution_clock::now();
	std::srand(time(0));

	//std::vector<double> x = { 1.0,0.01,0.1};
	//std::vector<double> lorenzParams = {10.0,8./3.0,28.0};
	//function f = function(function::lorenz_equation, lorenzParams, 3, 3);
	//std::vector<std::pair<std::vector<double>, double>> y = f.runge_kutta4th(make_pair(x, 0), 100, 0.001);
	//matrix::printMatrixToFile(y, "lorenz.txt");
	//plot3D("lorenz", "lorenz", "Sistema de Lorenz", "2:3:4");
	
	/////////////////////////////////////////////////////////////
	//// Dipolo magnético
	//std::vector<double> x = { 1.0,0.0 };
	//std::vector<double> pendulumParams = { 1.0,2.5};
	//
	//function magenticDipole = function(function::magenticDipole, pendulumParams, 2);
	//
	//std::vector<std::pair<std::vector<double>, double>> y = magenticDipole.runge_kutta4th(make_pair(x, 0), 100, 0.001);
	//matrix::printMatrixToFile(y, "./Dat/dipole_001_f2.5.txt");
	//
	//plotting::plot2D("dipole_001_f2.5", "./Dat/dipole_001_f2.5", "Magnetic dipole in an osciling magnectic field with frequency of 2.5", "1:2");
	///////////////////////////////////////////////////////////////

	std::vector<double> x = { 1.0,0.0 };
	std::vector<double> pendulumParams = { 1.0,0.0};
	double tf = 10.0 * PI;

	function magneticDipole = function(function::magenticDipole, make_pair(x, 0),pendulumParams, 2);
	
	std::vector<std::pair<std::vector<double>, double>> y = magneticDipole.runge_kutta4th(make_pair(x, 0), tf, 0.1);
	//matrix::printMatrixToFile(y, "./Dat/dipole_001_f2.5.txt");
	
	//plotting::plot2D("dipole_001_f2.5", "./Dat/dipole_001_f2.5", "Magnetic dipole in an osciling magnectic field with frequency of 2.5", "1:2");

	std::ofstream output;
	output.open("./Dat/TesteFourierNovo.txt");
	//for (unsigned int i = 0; i < 1e5; i++)
	//{
	//	fTransform = function::fourierTransform(y, i*0.01,0,tf);
	//	
	//	output << i*0.01 <<" " << std::norm(fTransform) << "\n";
	//}
	std::vector<std::pair<std::complex<double>, double>> fTransform = magneticDipole.fourierTransformRange(0.0, tf, 0.001, 0, 10);
	for (auto i:fTransform)
	{
		output << i.second  <<' ' << i.first.real() <<" " << i.first.imag() << " " << std::norm(i.first) << "\n";
	}
	output.close();
	std::vector<double> frequencies;

	for (unsigned int i = 0; i < fTransform.size(); i++)
	{
		if(std::norm(fTransform[i].first) > 0.8)
		{
			frequencies.emplace_back(std::norm(fTransform[i].first));
		}
	}
	std::cout << frequencies[round((frequencies.size()) / 2)] << "\n";
	//plotting::plot2D("fourierTransformDipole2_50", "./Dat/FourierTransform2_50", "Fourier Transform Dipole in a uniform magnetic field f 2.50", "1:2");


	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
	
	std::cout << "\n Time taken by the program: " << duration.count() << " miliseconds" << "\n";

	return 0;
}