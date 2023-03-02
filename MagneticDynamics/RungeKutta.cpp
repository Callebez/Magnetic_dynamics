#include "RungeKutta.hpp"
void RK4thSolver::runge_kutta_step(void (*func) (std::pair<double*, double>, double*, double*&),std::pair<double*, double> x,std::pair<double*, double>& y1 )
{
    // auto y1 = std::make_unique<std::vector<double>>(new double[rk_int->sysDim]);
    // y1->resize(rk_int->sysDim);
	// std::unique_ptr<double*> y1 = std::make_unique<double*> ();

	// *y1 = new double[rk_int->sysDim];
    // double* y1 = new double[rk_int->sysDim];
//    for (unsigned int j = 0; j < rk_int->sysDim; j++)
//             std::cout<<"x"<<j<<": " <<x.first[j]<<"\n";

    func(x,rk_int->params, k1);
    matrix::axpy(k1, x.first, rk_int->step / 2.0,rk_int->sysDim,k2aux);

    func(std::make_pair(k2aux, x.second + rk_int->step / 2.0), rk_int->params,k2);

    matrix::axpy(k2, x.first, rk_int->step / 2.0,rk_int->sysDim,k3aux);
    func(std::make_pair(k3aux, x.second + rk_int->step / 2.0), rk_int->params,k3);

    matrix::axpy(k3, x.first, rk_int->step ,rk_int->sysDim,k4aux);
    func(std::make_pair(k4aux, x.second + rk_int->step), rk_int->params, k4);
    y1.first = new double[rk_int->sysDim];
    for (unsigned int i = 0; i < rk_int->sysDim; i++)
    {
        y1.first[i] = x.first[i] + (rk_int->step /6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
    y1.second = x.second + rk_int->step;
    // for (unsigned int i = 0; i < rk_int->sysDim; i++)
        // std::cout<<y1.first[i]<<'\n';
} 
void RK4thSolver::RungeKuttaMethod(void (*func) (std::pair<double*, double>, double*, double*&), std::pair<double*, double> x0, double* param, double t_initial, double t_final, double h, uint dim)
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
    k1    = new double[rk_int->sysDim];
    k2    = new double[rk_int->sysDim]; 
    k3    = new double[rk_int->sysDim];
    k4    = new double[rk_int->sysDim];
    k2aux = new double[rk_int->sysDim]; 
    k3aux = new double[rk_int->sysDim];
    k4aux = new double[rk_int->sysDim];

    for (unsigned int i = 0; i < rk_int->n_iterations-1; i++)
    {
        runge_kutta_step(func, rk_int->data[i], rk_int->data[i+1]);
    }

    delete[] k1;
    delete[] k2;
    delete[] k3;
    delete[] k4;
    
    delete[] k2aux;
    delete[] k3aux;
    delete[] k4aux;
    // delete[] y1;

}
double* RK45thSolver::runge_kutta_step(void (*func) (std::pair<double*, double>, double*, double*&),std::pair<double*, double> x)
{
 
    erro.clear();
    erro.resize(rk_int->sysDim);
    y1 = new double[rk_int->sysDim];

    double* aux = new double[rk_int->sysDim];
    // std::cout<<k1[0]<<" "<<k1[1]<<" "<<k1[2]<<'\n';
    func(x,rk_int->params, k1);
    matrix::axpy(k1, x.first, rk_int->step * b21, rk_int->sysDim, aux);
    func(std::make_pair(aux, x.second + rk_int->step *a2), rk_int->params,k2);

    matrix::axpy(k1, x.first, rk_int->step * b31, rk_int->sysDim, k3aux);
    matrix::axpy(k2, k3aux  , rk_int->step * b32, rk_int->sysDim, aux);
    func(std::make_pair(aux, x.second + rk_int->step *a3), rk_int->params,k3);
    
    matrix::axpy(k1, x.first, rk_int->step * b41,rk_int->sysDim, k4aux1);
    matrix::axpy(k2, k4aux1,  rk_int->step * b42,rk_int->sysDim, k4aux2);
    matrix::axpy(k3, k4aux2 , rk_int->step * b43,rk_int->sysDim, aux);
    func(std::make_pair(aux, x.second + rk_int->step *a4), rk_int->params,k4);
    
    matrix::axpy(k1, x.first, rk_int->step * b51,rk_int->sysDim, k3aux);
    matrix::axpy(k2, k3aux,   rk_int->step * b52,rk_int->sysDim, k4aux1);
    matrix::axpy(k3, k4aux1,  rk_int->step * b53,rk_int->sysDim, k4aux2);
    matrix::axpy(k4, k4aux2,  rk_int->step * b54,rk_int->sysDim, aux);
    func(std::make_pair(aux, x.second + rk_int->step *a5), rk_int->params,k5);
    
    matrix::axpy(k1, x.first, rk_int->step * b61,rk_int->sysDim, k3aux);
    matrix::axpy(k2, k3aux,   rk_int->step * b62,rk_int->sysDim, k4aux1);
    matrix::axpy(k3, k4aux1,  rk_int->step * b63,rk_int->sysDim, k4aux2);
    matrix::axpy(k4, k4aux2,  rk_int->step * b64,rk_int->sysDim, k6aux);
    matrix::axpy(k5, k6aux,   rk_int->step * b65,rk_int->sysDim, aux);
    func(std::make_pair(aux, x.second + rk_int->step*a6), rk_int->params, k6);

    for (unsigned int i = 0; i < rk_int->sysDim; i++)
    {
        y1[i]   = x.first[i] + rk_int->step * (ch1 * k1[i] +  ch2 * k2[i] + ch3 * k3[i] + ch4 * k4[i] + ch5 * k5[i] + ch6 * k6[i] );
        erro[i] = ( rk_int->step *(d1*k1[i]+d3*k3[i] + d4*k4[i]+ d5*k5[i]+d6*k6[i])); 
    }

    delete[] aux;
    
    return y1;
} 
void RK45thSolver::RungeKuttaMethod(void (*func) (std::pair<double*, double>, double*, double*&), std::pair<double*, double> x0, double* param, double t_initial, double t_final, double h, uint dim)
{
    
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
    k1 = new double[rk_int->sysDim];
    k2 = new double[rk_int->sysDim]; 
    k3 = new double[rk_int->sysDim];
    k4 = new double[rk_int->sysDim];
    k5 = new double[rk_int->sysDim];
    k6 = new double[rk_int->sysDim];
     k3aux = new double[rk_int->sysDim];
    k4aux1= new double[rk_int->sysDim];
    k4aux2= new double[rk_int->sysDim];
    k6aux= new double[rk_int->sysDim];
    erro.clear();
    while (rk_int->data[i].second < t_final)
    {
        rk_int->data[i+1] = std::make_pair(runge_kutta_step(func, rk_int->data[i]), rk_int->data[i].second + rk_int->step);
        erroNorm = matrix::norm(erro);
        erroTotal+=erroNorm;
        h_new = (0.9*rk_int->step) * std::pow(((double)tol)/(double)(erroNorm), 1.0/5.0);

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

    delete[] k3aux;
    delete[] k4aux1;
    delete[] k4aux2;
    delete[] k6aux;
    delete[] k1;
    delete[] k2;
    delete[] k3;
    delete[] k4;
    delete[] k5;
    delete[] k6;
    delete[] y1;
}
