#include "RungeKutta.hpp"
double* rungeKutta::runge_kutta_step(double* (*func) (std::pair<double*, double>, double*),std::pair<double*, double> x)
{
    double* y1 = new double[rk_int->sysDim];
    double* k1 = new double[rk_int->sysDim];
    double* k2 = new double[rk_int->sysDim];
    double* k3 = new double[rk_int->sysDim];
    double* k4 = new double[rk_int->sysDim];

    k1 = func(x,rk_int->params);
    k2 = func(std::make_pair(matrix::axpy(k1, x.first, rk_int->step / 2.0,rk_int->sysDim), x.second + rk_int->step / 2.0), rk_int->params);
    k3 = func(std::make_pair(matrix::axpy(k2, x.first, rk_int->step / 2.0,rk_int->sysDim), x.second + rk_int->step / 2.0), rk_int->params);
    k4 = func(std::make_pair(matrix::axpy(k3, x.first, rk_int->step ,rk_int->sysDim), x.second + rk_int->step), rk_int->params);

    for (unsigned int i = 0; i < rk_int->sysDim; i++)
    {
        y1[i] = x.first[i] + (rk_int->step /6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
    return y1;
}
void rungeKutta::runge_kutta4th(double* (*func) (std::pair<double*, double>, double*), std::pair<double*, double> x0, double* param, double t_initial, double t_final, double h, uint dim)
{
    
    rk_int->n_iterations = (unsigned int)((t_final - t_initial) / h);
    
    rk_int->step = h;

    rk_int->t0 = t_initial;
    rk_int->tf = t_final;
    rk_int->params = param;
    rk_int->step = h; 
    rk_int->sysDim = dim;
    rk_int->data = std::vector<std::pair<double*,double>> (rk_int->n_iterations);
    rk_int->data[0] = x0; 

    for (unsigned int i = 0; i < rk_int->n_iterations-1; i++)
    {
        rk_int->data[i+1] = std::make_pair(runge_kutta_step(func, rk_int->data[i]), rk_int->data[i].second + h);
    }
}
