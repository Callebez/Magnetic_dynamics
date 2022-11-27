#include "functions.h"
function::function(std::vector<double>(*fun) (std::pair<std::vector<double>, double>, std::vector<double>), std::vector<double> param,unsigned int dim)
{
	n_dim = dim;
    params = param;
	func = fun;	
}

std::vector<double> function::lorenz_equation(std::pair<std::vector<double>, double>x, std::vector<double> param)
{
	std::vector<double> x_dot(3);
	x_dot[0] = param[0] * (x.first[1] - x.first[0]);
	x_dot[1] = x.first[0] * (param[2] - x.first[2]) - x.first[1];
	x_dot[2] = x.first[0] * x.first[1] - param[1] * x.first[2];
    
	return x_dot;
}
std::vector<double> function::sistemaAmigos(std::pair<std::vector<double>, double>x, std::vector<double> param)
{
    std::vector<double> x_dot(2);
    x_dot[0] = x.first[1];
    x_dot[1] = -param[0] + param[1]*x.first[0]*x.first[1];
    return x_dot;
}
std::vector<double> function::pendulum(std::pair<std::vector<double>, double>x, std::vector<double> param)
{
    std::vector<double> x_dot(2);
    x_dot[0] = x.first[1];
    x_dot[1] = -param[0] * sin(x.first[0]-sin(param[1]*x.second));
    return x_dot;
}
std::vector<double> function::runge_kutta_step(std::pair<std::vector<double>, double> x, double h)
{
    std::vector<double> y1(n_dim);
    std::vector<double> k1(n_dim);
    std::vector<double> k2(n_dim);
    std::vector<double> k3(n_dim);
    std::vector<double> k4(n_dim);

    k1 = func(x,params);
    k2 = func(std::make_pair(matrix::axpy(k1, x.first, h / 2.0), x.second + h / 2.0), params);
    k3 = func(std::make_pair(matrix::axpy(k2, x.first, h / 2.0), x.second + h / 2.0), params);
    k4 = func(std::make_pair(matrix::axpy(k3, x.first, h ), x.second + h), params);

    for (unsigned int i = 0; i < n_dim; i++)
    {
        y1[i] = x.first[i] + (h /6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
    return y1;
}
std::vector<std::pair<std::vector<double>, double>> function::runge_kutta4th(std::pair<std::vector<double>, double> x, double tf, double h)
{
    unsigned int n_iterations = (unsigned int)((tf - x.second) / h);
    std::vector<std::pair<std::vector<double>, double>> res(n_iterations);
    res[0] = x; 
    for (unsigned int i = 0; i < n_iterations-1; i++)
    {
        res[i+1] = std::make_pair(runge_kutta_step(res[i], h), res[i].second + h);
    }
    return res;
}
