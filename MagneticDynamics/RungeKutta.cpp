#include "RungeKutta.hpp"
double* RK4thSolver::runge_kutta_step(double* (*func) (std::pair<double*, double>, double*),std::pair<double*, double> x)
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
    free(k1);
    free(k2);
    free(k3);
    free(k4);
    return y1;
} 
void RK4thSolver::RungeKuttaMethod(double* (*func) (std::pair<double*, double>, double*), std::pair<double*, double> x0, double* param, double t_initial, double t_final, double h, uint dim)
{
    rk_int = std::make_unique<solution<double*>>();
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
double* RK45thSolver::runge_kutta_step(double* (*func) (std::pair<double*, double>, double*),std::pair<double*, double> x)
{
    double a2 = 0.25;
    double a3 = 0.375;
    double a4 = 12.0/13.0;
    double a5 = 1;
    double a6 = 0.5;
    double b21 = 0.25;
    double b31 = 3.0/32.0;
    double b32 = 9.0/32.0;
    double b41 = 1932.0/2197.0;
    double b42 = -7200.0/2197.0;
    double b43 = 7296.0/2197.0;
    double b51 = 439.0/216.0;
    double b52 = -8.0;
    double b53 = 3680.0/513.0;
    double b54 = -845.0/4104.0;
    double b61 = -8.0/27.0;
    double b62 = 2.0;
    double b63 = -3544.0/2565.0;
    double b64 = 1859.0/4104.0;
    double b65 = -11.0/40.0;
    // double c1 = 25.0/216.0;
    // double c3 = 1408.0/2565.0;
    // double c4 = 2197.0/4104.0;
    // double c5 = -0.20;
    double d1 = 1.0/360.0;
    double d3 = -128.0/4275.0;
    double d4 = -2197.0/75240.0;
    double d5 = 0.02;
    double d6 = 2.0/55.0;

    double ch1 = 16.0/135.0; 
    double ch2 = 0.0;
    double ch3 = 6656.0/12825.0	; 
    double ch4 = 28561.0/56430.0; 
    double ch5 = -9.0/50.0	; 
    double ch6 = 2.0/55.0; 
    erro.clear();
    erro.resize(rk_int->sysDim);
    // erro.clear()
    double* y1 = new double[rk_int->sysDim];
    double* k1 = new double[rk_int->sysDim];
    double* k2 = new double[rk_int->sysDim]; 
    double* k3 = new double[rk_int->sysDim];
    double* k4 = new double[rk_int->sysDim];
    double* k5 = new double[rk_int->sysDim];
    double* k6 = new double[rk_int->sysDim];
    
    double* k3aux= new double[rk_int->sysDim];
    double* k4aux1= new double[rk_int->sysDim];
    double* k4aux2= new double[rk_int->sysDim];
    double* k6aux= new double[rk_int->sysDim];

    k1 = func(x,rk_int->params);
    k2 = func(std::make_pair(matrix::axpy(k1, x.first, rk_int->step * b21, rk_int->sysDim), x.second + rk_int->step *a2), rk_int->params);
    k3aux = matrix::axpy(k1, x.first, rk_int->step * b31, rk_int->sysDim);
    k3 = func(std::make_pair(matrix::axpy(k2, k3aux, rk_int->step * b32, rk_int->sysDim), x.second + rk_int->step *a3), rk_int->params);
    k4aux1 = matrix::axpy(k1, x.first, rk_int->step * b41,rk_int->sysDim);
    k4aux2 = matrix::axpy(k2,k4aux1,rk_int->step*b42,rk_int->sysDim);
    k4 = func(std::make_pair(matrix::axpy(k3, k4aux2 ,rk_int->step * b43,rk_int->sysDim), x.second + rk_int->step*a4), rk_int->params);
    
    k3aux  = matrix::axpy(k1, x.first, rk_int->step * b51,rk_int->sysDim);
    k4aux1 = matrix::axpy(k2,k3aux,rk_int->step*b52,rk_int->sysDim);
    k4aux2 = matrix::axpy(k3,k4aux1, rk_int->step * b53,rk_int->sysDim);
    k5 = func(std::make_pair(matrix::axpy(k4, k4aux2, rk_int->step * b54,rk_int->sysDim), x.second + rk_int->step*a5), rk_int->params);
    
    k3aux  = matrix::axpy(k1, x.first, rk_int->step * b61,rk_int->sysDim);
    k4aux1 = matrix::axpy(k2,k3aux,rk_int->step*b62,rk_int->sysDim);
    k4aux2 = matrix::axpy(k3,k4aux1, rk_int->step * b63,rk_int->sysDim);
    k6aux = matrix::axpy(k4, k4aux2, rk_int->step * b64,rk_int->sysDim);
    k6 = func(std::make_pair(matrix::axpy(k5, k6aux,rk_int->step * b65,rk_int->sysDim), x.second + rk_int->step*a6), rk_int->params);
    
    for (unsigned int i = 0; i < rk_int->sysDim; i++)
    {
        y1[i]   = x.first[i] + rk_int->step * (ch1 * k1[i] +  ch2 * k2[i] + ch3 * k3[i] + ch4 * k4[i] + ch5 * k5[i] + ch6 * k6[i] );
        erro[i] = ( rk_int->step *(d1*k1[i]+d3*k3[i] + d4*k4[i]+ d5*k5[i]+d6*k6[i])); 
    }

    free(k1);
    free(k2);
    free(k3);
    free(k3aux);
    free(k4);
    free(k4aux1);
    free(k4aux2);
    free(k5);
    free(k6aux);
    free(k6);
    return y1;
} 
void RK45thSolver::RungeKuttaMethod(double* (*func) (std::pair<double*, double>, double*), std::pair<double*, double> x0, double* param, double t_initial, double t_final, double h, uint dim)
{
    erro.clear();
    rk_int = std::make_unique<solution<double*>>();
    rk_int->n_iterations = (unsigned int)((t_final - t_initial) / h);  
    rk_int->step = h;
    rk_int->t0 = t_initial;
    rk_int->tf = t_final;      
    rk_int->params = param;
    rk_int->step = h; 
    rk_int->sysDim = dim;
    rk_int->data = std::vector<std::pair<double*,double>>(2);
    rk_int->data[0] = x0; 
    double erroNorm =0;
    double h_new = h;
    double tol = 1e-8;
    double erroTotal = 0;
    uint i = 0;
 
    while (rk_int->data[i].second < t_final)
    {
        rk_int->data[i+1] = std::make_pair(runge_kutta_step(func, rk_int->data[i]), rk_int->data[i].second + rk_int->step);
        erroNorm = matrix::norm(erro);
        erroTotal+=erroNorm;
        h_new = rk_int->step*std::pow(((double)tol)/(double)(erroNorm), 1.0/5.0);

        while(erroNorm>=tol)
        {
            rk_int->step=h_new;
            rk_int->data[i+1] = std::make_pair(runge_kutta_step(func, rk_int->data[i]), rk_int->data[i].second + rk_int->step);
            h_new = rk_int->step*std::pow(((double)tol)/(double)(erroNorm), 1.0/5.0);
            erroNorm = matrix::norm(erro);        
        }
        i++;
        rk_int->data.resize(i+2);   
    }
}
