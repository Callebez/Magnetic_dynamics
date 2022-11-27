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
	std::vector<double> pendulumParams = { 1.0,0.75};
	
	function magenticDipole = function(function::magenticDipole, pendulumParams, 2);
	
	std::vector<std::pair<std::vector<double>, double>> y = magenticDipole.runge_kutta4th(make_pair(x, 0), 100.0 * PI, 0.1);
	matrix::printMatrixToFile(y, "./Dat/dipole_002_f0.75.txt");
	
	plotting::plot2D("dipole_002_f0.75", "./Dat/dipole_002_f0.75", "Magnetic dipole in an osciling magnectic field with frequency of 0.75", "1:2");

	std::vector <std::complex<double>> yComplex;
	std::complex<double> integral;
	std::ofstream output;

	output.open("./Dat/FourierTransform0_75.txt");
	for (unsigned int i = 0; i < 1e3; i++)
	{
		yComplex = function::fourierTransform(y, i*0.1);
		integral = function::simpsonOneThird(yComplex, 0, 100.0 * PI);
		output << i*0.1 <<" " << std::norm(integral) << "\n";
	}
	output.close();
	plotting::plot2D("fourierTransformDipole0_75", "./Dat/FourierTransform0_75", "Fourier Transform Dipole in a uniform magnetic field f 0.75", "1:2");


	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
	
	std::cout << "\n Time taken by the program: " << duration.count() << " miliseconds" << "\n";

	return 0;
}