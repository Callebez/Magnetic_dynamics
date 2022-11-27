#include <iostream>
#include<fstream>
#include<chrono>
#include<cstdlib>
#include "matrices.h"
#include "functions.h"
#include <ctime>
#include "plotting.h"


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

	std::vector<double> x = { 1.0,0.0 };
	std::vector<double> pendulumParams = { 1.0,2.5};
	function pendulum = function(function::pendulum, pendulumParams, 2);
	std::vector<std::pair<std::vector<double>, double>> y = pendulum.runge_kutta4th(make_pair(x, 0), 100, 0.001);
	matrix::printMatrixToFile(y, "./Dat/dipole_001_f2.5.txt");

	plotting::plot2D("dipole_001_f2.5", "./Dat/dipole_001_f2.5", "Magnetic dipole in an osciling magnectic field with frequency of 2.5", "1:2");
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
	std::cout << "\n Time taken by function: " << duration.count() << " miliseconds" << "\n";

	return 0;
}