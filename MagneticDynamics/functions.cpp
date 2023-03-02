#include "functions.hpp"

void MagneticDipole::magenticDipole(std::pair<double*, double> x, double* param, double*& res)
{   
    res[0] = x.first[1];
    // angle - twoPi * floor( angle / twoPi )
    res[1] = -param[0] * sin((x.first[0] - 2*3.14159265359*floor(x.first[0]/2*3.14159265359)) - 3.14159265359/2.0*sin(param[1]*x.second));
}
void Lorenz::lorenz_equation(std::pair<double*, double>x, double* param, double*& res)
{
    res[0] = param[0] * (x.first[1] - x.first[0]);   
    res[1] = x.first[0] * (param[2] - x.first[2]) - x.first[1]; 
    res[2] = x.first[0] * x.first[1] - param[1] * x.first[2];    
};
void Function::applyFourierTransform(double initialFrequency, double finalFrequency, double frequencyStep)
{
    std::unique_ptr<fourier> F = std::make_unique<fourier>();   
    F->fourierFrequencySpectrumAbsoluteValue(solver->rk_int,initialFrequency,finalFrequency,frequencyStep);
    fourier_transform = F->get_FrequencySpectrum();

}

// solution<double*>* Function::createSignal(function signalFunc, double* params_, double t0_, double tf_, double step_, uint sysDim_)
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

void TripleMagneticDipole::tripleMagenticDipole(std::pair<double *, double> x, double *param, double *&res)
{
    double M1 = 1.0;
    double M2 = 1.0;
    double M3 = 1.0;
    double consts = 1.0;
    res[0] = x.first[3];
    res[1] = x.first[4];
    res[2] = x.first[5];
    res[3] = consts*param[0]*M1*(3.0*sqrt(3.0)*M2*cos(x.first[0]+x.first[1])-3.0*sqrt(3.0)*M3*cos(x.first[0]+x.first[2]) + M2*cos(x.first[1])*sin(x.first[0])+M3*cos(x.first[2])*sin(x.first[0])+5.0*M2*cos(x.first[0])*sin(x.first[0])+ 5.0*M3*cos(x.first[0])*sin(+x.first[2]));
    res[4] = consts*param[1]*M2*(3.0*sqrt(3.0)*M1*cos(x.first[0]+x.first[1])+ (M1*cos(x.first[0])-8.0*M3*cos(x.first[2]))*sin(x.first[1])+cos(x.first[1])*(5.0*M1*sin(x.first[0]) - 4.0*M3*sin(x.first[3])));
    res[5] = -consts*param[2]*M3*(3.0*sqrt(3.0)*M1*cos(x.first[0]+x.first[2])+ cos(x.first[2])*(-5.0*M1*sin(x.first[0]) + 4.0*M2*sin(x.first[1]))+sin(x.first[2])*(-M1*cos(x.first[0])+8.0*M2*cos(x.first[1])));
}
