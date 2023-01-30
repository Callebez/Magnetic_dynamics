#include "functions.hpp"


// void function::applyRunge_kutta4th(std::pair<double*, double> x0, double t_initial, double t_final, double h)
// {
//     std::unique_ptr<rungeKutta> Rk4 = std::make_unique<rungeKutta>();
//     // rungeKutta* Rk4 = new rungeKutta;
//     // rk_int = std::make_unique<solution<double*>>;
//     // rk_int  = new solution<double*>;
//     Rk4->runge_kutta4th(func,x0,params,t_initial,t_final,h, x_dim);
//     rk_int = Rk4->getSolution();
//     // Rk4->rk_int->data.clear();
//     // delete Rk4; 
//     // return *rk_int;
// } 
// double* function::magenticDipole(std::pair<double*, double> x, double* param)
// {   
//     double* x_dot = new double[2];
//     x_dot[0] = x.first[1];
//     // x_dot[1] = -param[0] * param[1]*x.first[0];
//     x_dot[1] = -param[0] * sin(x.first[0]-sin(param[1]*x.second));
//     return x_dot;
// }
void Lorenz::lorenz_equation(std::pair<double*, double>x, double* param, double*& res)
{
    // std::cout<<"x : "<< x.first[0]<<" "<<x.first[1]<<" "<< x.first[2]<<"\n";
    res[0] = param[0] * (x.first[1] - x.first[0]);
    res[1] = x.first[0] * (param[2] - x.first[2]) - x.first[1];
    res[2] = x.first[0] * x.first[1] - param[1] * x.first[2];    
    // std::cout<<param[0] * (x.first[1] - x.first[0])<<"\n";//<< '  res: '<< res[0]<<"\n, ";
    // std::cout<<x.first[0] * (param[2] - x.first[2]) - x.first[1]<<"\n res: " << res[1]<<"\n , ";
    // std::cout<<x.first[0] * x.first[1] - param[1] * x.first[2]<<"res :" <<res[2]<<"\n";
};
// void Functions::applySolver(Solver& solverMethod, std::pair<double*,double> x0,double t0, double tf, double step)
// {

//     ClientCode(solverMethod, func, x0, t0, tf, step);
// }

// void function::applyFourierTransform(double initialFrequency, double finalFrequency, double frequencyStep)
// {
//     std::unique_ptr<fourier> F = std::make_unique<fourier>();   
//     // auto fourier_transform = std::make_unique<solution<double>>;
//     F->fourierFrequencySpectrumAbsoluteValue(rk_int,initialFrequency,finalFrequency,frequencyStep);
//     fourier_transform = F->get_FrequencySpectrum();
//     // delete F;
//     // return fourier_transform;
// }

// solution<double*>* function::createSignal(function signalFunc, double* params_, double t0_, double tf_, double step_, uint sysDim_)
// {
// 	solution<double*>* signal = new solution<double*>;

// 	uint n = (uint)(tf_- t0_)/step_;
// 	signal->data = std::vector<std::pair<double*,double>> (n);
	
// 	double* x = new double[sysDim_];

// 	for(uint i = 0; i < n; i++)
// 	{
// 		signal->data[i].first = new double[sysDim_];
// 	}

// 	signal->t0 = t0_;
// 	signal->tf = tf_;
// 	signal->n_iterations = n;
// 	signal->sysDim = sysDim_;
	
// 	for(uint i = 0; i < n; i++)
// 	{
// 		signal->data[i].first = signalFunc.func(std::make_pair(x,t0_ + i*step_), params_);
// 		signal->data[i].second = t0_ + i*step_;
// 	}	

//     return signal;
// }
