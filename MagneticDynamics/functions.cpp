#include "functions.hpp"


solution<double*>* function::applyRunge_kutta4th(std::pair<double*, double> x0, double t_initial, double t_final, double h)
{
    rungeKutta Rk4;

    Rk4.runge_kutta4th(func,x0,params,t_initial,t_final,h, x_dim);
    rk_int = Rk4.getSolution();
    return rk_int;
}
double* function::magenticDipole(std::pair<double*, double> x, double* param)
{   
    double* x_dot = new double[2];
    x_dot[0] = x.first[1];
    // x_dot[1] = -param[0] * param[1]*x.first[0];
    x_dot[1] = -param[0] * sin(x.first[0]-sin(param[1]*x.second));
    return x_dot;
}
double* function::lorenz_equation(std::pair<double*, double>x, double* param)
{
    double* x_dot = new double[3];
    x_dot[0] = param[0] * (x.first[1] - x.first[0]);
    x_dot[1] = x.first[0] * (param[2] - x.first[2]) - x.first[1];
    x_dot[2] = x.first[0] * x.first[1] - param[1] * x.first[2];
    
    return x_dot;
};

solution<double>* function::applyFourierTransform(double initialFrequency, double finalFrequency, double frequencyStep)
{
    // std::cout<<"memory of the rk4 "<<&rk_int<<"\n";
    fourier F;
    // std::cout<<"memory of the function "<<&(F.fourierTransform)<<"\n";
    
    F.fourierFrequencySpectrumAbsoluteValue(rk_int,initialFrequency,finalFrequency,frequencyStep);
    fourier_transform = F.get_FrequencySpectrum();
    return fourier_transform;
}

solution<double*>* function::createSignal(function signalFunc, double* params_, double t0_, double tf_, double step_, uint sysDim_)
{
	solution<double*>* signal = new solution<double*>;

	uint n = (uint)(tf_- t0_)/step_;
	signal->data = new std::pair<double*,double>[n];
	
	double* x = new double[sysDim_];

	for(uint i = 0; i < n; i++)
	{
		signal->data[i].first = new double[sysDim_];
	}

	signal->t0 = t0_;
	signal->tf = tf_;
	signal->n_iterations = n;
	signal->sysDim = sysDim_;
	
	for(uint i = 0; i < n; i++)
	{
		signal->data[i].first = signalFunc.func(std::make_pair(x,t0_ + i*step_), params_);
		signal->data[i].second = t0_ + i*step_;
	}	

    return signal;
}



















// function::function(double* (*fun) (std::pair<double*, double>, double*), std::pair<double*, double> x0, double* param, unsigned int dim)
// {
// 	n_dim = dim;
//     params = param;
// 	func = fun;	
// }

// double* lorenz::lorenz_equation(std::pair<double*, double>x, std::vector<double> param)
// {
// 	double* x_dot = new double[3];
// 	x_dot[0] = param[0] * (x.first[1] - x.first[0]);
// 	x_dot[1] = x.first[0] * (param[2] - x.first[2]) - x.first[1];
// 	x_dot[2] = x.first[0] * x.first[1] - param[1] * x.first[2];
    
// 	return x_dot;
// }
// double* function::sistemaAmigos(std::pair<double*, double>x, std::vector<double> param)
// {
//     double* x_dot = new double[2];
//     x_dot[0] = x.first[1];
//     x_dot[1] = -param[0] + param[1]*x.first[0]*x.first[1];
//     return x_dot;
// }
// double* function::magenticDipole( std::pair<double*, double> x, double* param)
// {   
//     double* x_dot = new double[2];
//     x_dot[0] = x.first[1];
//     // x_dot[1] = -param[0] * param[1]*x.first[0];
//     x_dot[1] = -param[0] * sin(x.first[0]-sin(param[1]*x.second));
//     return x_dot;
// }
// // double* function::runge_kutta_step(std::pair<double*, double> x, double h)
// // {
// //     double* y1 = new double[n_dim];
// //     double* k1 = new double[n_dim];
// //     double* k2 = new double[n_dim];
// //     double* k3 = new double[n_dim];
// //     double* k4 = new double[n_dim];

// //     k1 = func(x,params);
// //     k2 = func(std::make_pair(matrix::axpy(k1, x.first, h / 2.0,n_dim), x.second + h / 2.0), params);
// //     k3 = func(std::make_pair(matrix::axpy(k2, x.first, h / 2.0,n_dim), x.second + h / 2.0), params);
// //     k4 = func(std::make_pair(matrix::axpy(k3, x.first, h ,n_dim), x.second + h), params);

// //     for (unsigned int i = 0; i < n_dim; i++)
// //     {
// //         y1[i] = x.first[i] + (h /6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
// //     }
// //     return y1;
// // }
// // solution* function::runge_kutta4th(std::pair<double*, double> x, double tf, double h)
// // {
// //     unsigned int n_iterations = (unsigned int)((tf - x.second) / h);
// //     std::pair<double*,double>* res= new std::pair<double*,double>[n_iterations];
// //     res[0] = x; 

// //     for (unsigned int i = 0; i < n_iterations-1; i++)
// //     {
// //         res[i+1] = std::make_pair(runge_kutta_step(res[i], h), res[i].second + h);
// //     }

// //     solved->data = res;
// //     solved->size = n_iterations;
// //     // solved->x0 = x;
// //     solved->step = h;
// //     solved->sysDim = n_dim;
// //     // free(res);
// //     return solved;
// // }

// std::complex<double> complexExp(double w, double t)
// {
//     return std::complex<double>(cos(w * t), -sin(w * t));
// }
// std::complex<double> multiRealByComplex(std::complex<double> z, double x)
// {
//     return std::complex<double>(z.real() * x, z.imag() * x);
// }
// std::complex<double> function::fourierTransform(double frequency)
// {
//     double t0 = solved->data[0].second;
//     double tf = solved->data[solved->n_iterations-1].second;
//     std::complex<double>* fexp = new std::complex<double> [solved->n_iterations];
//     std::complex<double> transform;
//     for (unsigned int i = 0; i < solved->n_iterations; i ++)
//     {
//         fexp[i] = complexExp(frequency, solved->data[i].second)* solved->data[i].first[0];
//     }
//     transform = 2.0/(tf-t0)*simpsonOneThird(fexp, t0, tf ,solved->n_iterations);
//     // std::cout<<params[0]<<" "<<params[1]<<"\n";
//     return transform;
    
// }

// std::complex<double> function::simpsonOneThird(std::complex<double>* f,double t0,double tf, uint size)
// {
//     unsigned int n;
//     if (size % 2 != 0) { n = size - 1; }
//     else { n = size;}
//     std::complex<double> h = (tf - t0) / (double)size;
//     std::complex<double> sum = {0.0, 0.0};
//     sum += f[0]*h/3.0;
//     sum += f[size]* h / 3.0;
//     for (unsigned int i = 1; i < uint((n-1)/2); i++)
//     {
//         sum += f[2 * i - 1] * 4.0 * h / 3.0;
//         sum += f[2 * i]* 2.0 * h / 3.0;
//     }
//     return sum;
// }
// void function::printRungeKutta(char* filename)
// {
//     std::fstream plot; 
// 	plot.open(filename,std::fstream::out);
// 	for(uint i = 0; i < solved->n_iterations; i++)
// 	{
// 		plot<<solved->data[i].second;
// 		for(uint j = 0; j < solved->sysDim; j++)
// 		{
// 			plot<<" "<<solved->data[i].first[j];
// 		}
// 		plot<<"\n";
// 	}
// 	plot.close();
// }

// fourier function:: fourierTransformRange(double step, double f0, double ff)
// {
//     uint n_points = uint((ff-f0)/step);
//     fourierT->size = n_points;
//     fourierT->step = step;
//     fourierT->data.resize(n_points);
//     double f_aux;

// 	for(int i = 0 ; i < n_points; i++)
// 	{
//         f_aux =  std::pow(std::norm(fourierTransform(i*step)),2);

// 		fourierT->data[i].first = f_aux;
// 		fourierT->data[i].second = i*step;
// 		// frequencies<<fourier->data[i].first << " " <<fourier->data[i].second <<std::endl;
// 	}
//     return *fourierT;
// }
// void function::printFourierTransformRange(char* filename)
// {
//     std::fstream plot; 
// 	plot.open(filename,std::fstream::out);
//     for(uint i = 0; i < fourierT->size; i++)
//     {
//         plot<<fourierT->data[i].second <<" "<<fourierT->data[i].first <<std::endl; 
//     }
//     plot.close();
// }
// std::vector<std::pair<double,double>> function::findFrquenciesFourier(double tresholdAmplitude)
// {
//     std::vector<std::pair<double,double>> frequencies;
//     std::vector<std::pair<double,double>> frequenciesAux;
    
//     for(int i = 0; i < fourierT->size; i++)
// 	{
// 		while(fourierT->data[i].first < tresholdAmplitude && i < fourierT->size)
// 		{ i++;}
// 		while (fourierT->data[i].first >= tresholdAmplitude )
// 		{
// 			frequencies.emplace_back(fourierT->data[i]);
// 			i++;
//             // std::cout<<frequencies.back().first<<" "<<'frequencies.back().second<<"\n";
// 		}
// 		if(frequencies.size() > 1)
// 		{
// 			frequenciesAux.emplace_back(peakSearch(frequencies,0,frequencies.size()));
// 			frequencies.clear();
// 		}
// 	}
//     return frequenciesAux;
// }
// std::pair<double,double> function::peakSearch(std::vector<std::pair<double,double>> v, int start, int end)
// {
// 	int mid; 
// 	mid = (end+start+1)/2;
// 	if(((v[mid].first > v[mid+1].first)&& mid == start)|| (v[mid].first > v[mid-1].first && mid == end))
// 	{
// 		return v[mid];
// 	}
// 	else if(v[mid].first > v[mid-1].first && v[mid].first > v[mid+1].first)
// 	{
// 		return v[mid];
// 	}
// 	// If right neighbor is higher then right subpart must have a peak.
// 	else if(v[mid].first <= v[mid+1].first)
// 	{
// 		return peakSearch(v, mid+1, end);
// 	}
// 	// If left neighbor is higher then left subpart must have a peak.
// 	else if(v[mid].first <= v[mid-1].first)
// 	{
// 		return peakSearch(v, start,mid-1);
// 	}
// 	return std::make_pair(0.0,0.0);
// }
