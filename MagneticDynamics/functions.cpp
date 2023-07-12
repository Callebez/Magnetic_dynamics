#include "functions.hpp"
double atr = 0.0;
void MagneticDipole::magenticDipole(std::pair<double*, double> x, double* param, double*& res)
{   
    res[0] = x.first[1];
    res[1] = -sin(x.first[0] -param[0]*sin(param[1]*x.second))- atr*x.first[1]; //


    //First approximation 
    // res[0] = x.first[1]; 
    // res[1] = -x.first[0] + param[0]*std::pow(param[1],2)*sin(param[1]*x.second); 
    
    // Duffing approximation 
    // res[0] = x.first[1]; 
    // res[1] = -x.first[0] + std::pow(x.first[0],3)/6.0 + param[0]*std::pow(param[1],2)*sin(param[1]*x.second); 
    

    //Using Hamilton Equations
    // res[0] = x.first[1] - param[0]*param[1]*cos(param[1]*x.second);
    // res[1] =-1*sin(x.first[0]);
}
void MagneticDipoleHarmonic::magneticDipoleHarmonic(std::pair<double*, double> x, double* param, double*& res)
{
    double sol = param[0]*std::pow(param[1],2)*sin(param[1]*x.second)/(std::pow(param[1],2)-1.0) + cos(x.second)+std::pow(param[1],3)*param[0]/(param[1]*param[1]-1.0)*sin(x.second);
    res[0] = x.first[1]; 
    res[1] = -x.first[0] + param[0]*std::pow(param[1],2)*sin(param[1]*x.second) + sol- sin(sol);    
}
void MagneticDipoleDuffing::magneticDipoleDuffing(std::pair<double*, double> x, double* param, double*& res)
{
    res[0] = x.first[1]; 
    res[1] = -x.first[0] + param[0]*sin(param[1]*x.second)+ std::pow(x.first[0] - param[0]*sin(param[1]*x.second),3)/6.0 + param[0]*std::pow(param[1],2)*sin(param[1]*x.second)- atr*x.first[1]; 
}
void Lorenz::lorenz_equation(std::pair<double*, double>x, double* param, double*& res)
{   
    res[0] = param[0]*(x.first[1]-x.first[0]);
    res[1] = x.first[0]*(param[1]-x.first[2])-x.first[1];
    res[2] = x.first[0]*x.first[1] - x.first[2]*param[2];
};
void Function::applyFourierTransform(double initialFrequency, double finalFrequency, double frequencyStep)
{
    std::unique_ptr<fourier> F = std::make_unique<fourier>();   
    F->fourierFrequencySpectrumAbsoluteValue(solver->rk_int,initialFrequency,finalFrequency,frequencyStep);
    fourier_transform = F->get_FrequencySpectrum();

}

// solution<double*>* Function::createSignal(Function signalFunc, double* params_, double t0_, double tf_, double step_, uint sysDim_)
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
//         signalFunc.func(std::make_pair(x,t0_ + i*step_), params_, x);
// 		signal->data[i].first = x;
// 		signal->data[i].second = t0_ + i*step_;
// 	}	

//     return signal;
// }

void TripleMagneticDipole::tripleMagenticDipole(std::pair<double *, double> x, double *param, double *&res)
{
    double M1 = param[0];
    double M2 = param[1];
    double M3 = param[2];
    double consts = 1.0;
    res[0] = x.first[3];
    res[1] = x.first[4];
    res[2] = x.first[5];
    res[3] = -(consts*M1*(3.0*sqrt(3.0)*M2*cos(x.first[0] + x.first[1]) - 3.0*sqrt(3.0)*M3*cos(x.first[0] + x.first[2]) + 5.0*M2*cos(x.first[1])*sin(x.first[0]) + 5.0*M3*cos(x.first[2])*sin(x.first[0]) + M2*cos(x.first[0])*sin(x.first[1]) + M3*cos(x.first[0])*sin(x.first[2])))/(36.*sqrt(3.0));
    res[4] = -(consts*M2*(3.0*sqrt(3.0)*M1*cos(x.first[0] + x.first[1]) + (5.0*M1*cos(x.first[0]) - 4.0*M3*cos(x.first[2]))*sin(x.first[1]) + cos(x.first[1])*(M1*sin(x.first[0]) - 8.0*M3*sin(x.first[2]))))/(36.*sqrt(3.0));
    res[5] = (consts*M3*(3.0*sqrt(3.0)*M1*cos(x.first[0] + x.first[2]) + cos(x.first[2])*(-(M1*sin(x.first[0])) + 8.0*M2*sin(x.first[1])) + (-5.0*M1*cos(x.first[0]) + 4.0*M2*cos(x.first[1]))*sin(x.first[2])))/(36.*sqrt(3.0));
}
void TripleMagneticDipole::tripleMagenticDipoleAtr(std::pair<double *, double> x, double *param, double *&res)
{
    double M1 = param[0];
    double M2 = param[1];
    double M3 = param[2];
    double consts = 1.0;
    double atr = 0.05;
    res[0] = x.first[3];
    res[1] = x.first[4];
    res[2] = x.first[5];
    res[3] = -(consts*M1*(3.0*sqrt(3.0)*M2*cos(x.first[0] + x.first[1]) - 3.0*sqrt(3.0)*M3*cos(x.first[0] + x.first[2]) + 5.0*M2*cos(x.first[1])*sin(x.first[0]) + 5.0*M3*cos(x.first[2])*sin(x.first[0]) + M2*cos(x.first[0])*sin(x.first[1]) + M3*cos(x.first[0])*sin(x.first[2])))/(36.*sqrt(3.0)) - atr*x.first[3];
    res[4] = -(consts*M2*(3.0*sqrt(3.0)*M1*cos(x.first[0] + x.first[1]) + (5.0*M1*cos(x.first[0]) - 4.0*M3*cos(x.first[2]))*sin(x.first[1]) + cos(x.first[1])*(M1*sin(x.first[0]) - 8.0*M3*sin(x.first[2]))))/(36.*sqrt(3.0))- atr*x.first[4];
    res[5] = (consts*M3*(3.0*sqrt(3.0)*M1*cos(x.first[0] + x.first[2]) + cos(x.first[2])*(-(M1*sin(x.first[0])) + 8.0*M2*sin(x.first[1])) + (-5.0*M1*cos(x.first[0]) + 4.0*M2*cos(x.first[1]))*sin(x.first[2])))/(36.*sqrt(3.0))- atr*x.first[5];
}
