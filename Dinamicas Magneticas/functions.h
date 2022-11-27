#pragma once
#include "matrices.h"
class function
{
	private:
		
		unsigned int n_dim;
		std::vector<double> params;
		std::vector<double> (*func) (std::pair<std::vector<double>, double>, std::vector<double>);// F(x,t,param) 

	public:
		std::vector<double> operator ()(std::pair<std::vector<double>, double> x)
		{
			return func(x, params);
		};
		function(std::vector<double>(*fun) (std::pair<std::vector<double>, double>, std::vector<double>), std::vector<double> param, unsigned int n_dim);
	
		std::vector<double> runge_kutta_step(std::pair<std::vector<double>, double> x, double h);
		std::vector<std::pair<std::vector<double>, double>> runge_kutta4th(std::pair<std::vector<double>, double> x, double tf, double h);

		static std::vector<double> lorenz_equation(std::pair<std::vector<double>, double>x, std::vector<double> param);
		static std::vector<double> magenticDipole(std::pair<std::vector<double>, double>x, std::vector<double> param);
		static std::vector<double> sistemaAmigos(std::pair<std::vector<double>, double>x, std::vector<double> param);
		static std::complex<double> simpsonOneThird(std::vector<std::complex<double>>f, double t0, double tf);
		static std::vector <std::complex<double>> fourierTransform(std::vector<std::pair<std::vector<double>, double>> f, double w);

		static std::vector<std::pair<std::vector<double>, double>> func_test(unsigned n_points, double t0, double tf);


};

