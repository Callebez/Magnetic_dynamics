// #pragma once
// #include <complex>
// #include<cstdlib>

// #include "functions.hpp"
// #include "Solutions.hpp"

// class integration
// {
// private: 
//     std::complex<double> integral;
// public:
// /**
//  * @brief Function to evaluate the integral of a given function.
//  * 
//  * Evaluate the integral using the Simpson's one third rule (https://en.wikipedia.org/wiki/Simpson%27s_rule). It takes not a function as argument, but a array of complex numbers, 
//  * that is because this function it is only used in conjunction with the Fourier Transform.  
//  * 
//  * @param f Array of complex numbers storing the evaluation of a function; 
//  * @param t0 Lower limit of integration;
//  * @param tf Upper limite of integratino; 
//  * @param size Size of the array of complex numbers. 
//  * @return std::complex<double> Value of the integral
//  */
//     static std::complex<double> simpsonOneThird(std::vector<std::complex<double>> f,double t0,double tf, uint size);
//     inline std::complex<double> get_integral(){return integral;};
// };
