#pragma once
#include "matrices.h"
class function
{
	private:
		
		unsigned int x_dim;
		unsigned int f_dim;
		std::vector<double> params;
		std::vector<double> (*func) (std::pair<std::vector<double>, double>, std::vector<double>);// F(x,t,param) 

	public:
		std::vector<double> operator ()(std::pair<std::vector<double>, double> x)
		{
			return func(x, params);
		};
		function(std::vector<double>(*fun) (std::pair<std::vector<double>, double>, std::vector<double>), std::vector<double> param, unsigned dom_dim, unsigned imag_dim);

		static std::vector<double> runge_kutta_step(function f, std::pair<std::vector<double>, double> x, double h);
		static std::vector<std::pair<std::vector<double>, double>> runge_kutta4th(function f, std::pair<std::vector<double>, double> x, double tf, double h);

		static std::vector<double> lorenz_equation(std::pair<std::vector<double>, double>x, std::vector<double> param);
		static std::vector<double> pendulum(std::pair<std::vector<double>, double>x, std::vector<double> param);
		static std::vector<double> sistemaAmigos(std::pair<std::vector<double>, double>x, std::vector<double> param);

};

