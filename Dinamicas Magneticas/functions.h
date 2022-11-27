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
		static std::vector<double> pendulum(std::pair<std::vector<double>, double>x, std::vector<double> param);
		static std::vector<double> sistemaAmigos(std::pair<std::vector<double>, double>x, std::vector<double> param);
		//std::vector<long double> lyapunov(std::vector<std::vector<double>>(*jacobian)(std::pair<std::vector<double>, double>), std::pair<std::vector<double>, double> x, double tf, double step);


};

