#include "functions.h"
function::function(std::vector<double>(*fun) (std::pair<std::vector<double>, double>, std::vector<double>), std::pair<std::vector<double>, double> x0, std::vector<double> param,unsigned int dim)
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
std::vector<double> function::magenticDipole(std::pair<std::vector<double>, double>x, std::vector<double> param)
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
    solved = res;
    x0 = x;
    return res;
}
std::vector<std::pair<std::vector<double>, double>> function::func_test(unsigned n_points, double t0, double tf)
{
    std::vector<std::pair<std::vector<double>, double>> f(n_points+1);
    double t = t0;
    double h = (tf-t0)/n_points;
    std::vector<double> x(1);
    for (unsigned int i = 0; i < n_points+1; i++)
    {

        x = {sin(2.0*t)};
        f[i] = std::make_pair(x, t); 
        t += h;
        x.clear();
    }
    return f;
}
std::complex<double> complexExp(double w, double t)
{
    return std::complex<double>(cos(w * t), -sin(w * t));
}
std::complex<double> multiRealByComplex(std::complex<double> z, double x)
{
    return std::complex<double>(z.real() * x, z.imag() * x);
}
std::complex<double> function::fourierTransform(std::vector<std::pair<std::vector<double>, double>> f,double frequency, double t0, double tf)
{
    unsigned int n = f.size();
    std::vector<std::complex<double>> fexp(n);
    std::complex<double> transform;
    for (unsigned int i = 0; i < n; i ++)
    {
        fexp[i] = complexExp(frequency, f[i].second)* f[i].first[0];
    }
    transform = 2.0/(tf-t0)*simpsonOneThird(fexp, t0, tf);
    return transform;
}
std::vector<std::complex<double>> function::fourierTransformRange(std::vector<std::pair<std::vector<double>, double>> f, double inicialFrequency, double finalFrequency, double step, double t0, double tf)
{
    unsigned int n_iterations = (unsigned int((finalFrequency - inicialFrequency) / step));
    std::vector<std::complex<double>> fTransform (n_iterations);
    for (unsigned int i = 0; i < n_iterations; i++)
    {
        fTransform[i] = function::fourierTransform(f, i * step, 0, tf);
    }
    return fTransform;
}
std::vector<std::pair<std::complex<double>, double>> function::fourierTransformRange(double inicialFrequency, double finalFrequency, double step, double t0, double tf)
{
    unsigned int n_iterations = (unsigned int((finalFrequency - inicialFrequency) / step));
    std::vector<std::pair<std::complex<double>, double>> fTransform(n_iterations);
    for (unsigned int i = 0; i < n_iterations; i++)
    {
        fTransform[i].first = function::fourierTransform(solved, i * step, 0, tf);
        fTransform[i].second = i * step;

    }
    return fTransform;
}
std::complex<double> function::simpsonOneThird(std::vector<std::complex<double>> f,double t0,double tf)
{
    unsigned int n;
    if (f.size() % 2 != 0) { n = f.size() - 1; }
    else { n = f.size();}
    std::complex<double> h = (tf - t0) / (double)f.size();
    std::complex<double> sum = {0.0, 0.0};
    sum += f[0]*h/3.0;
    sum += f.back()* h / 3.0;
    for (unsigned int i = 1; i < int((n-1)/2); i++)
    {
        sum += f[2 * i - 1] * 4.0 * h / 3.0;
        sum += f[2 * i]* 2.0 * h / 3.0;
    }
    return sum;
}
