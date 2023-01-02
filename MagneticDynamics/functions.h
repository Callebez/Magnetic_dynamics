// #pragma once
// #include "matrices.h"
// #include "RungeKutta.hpp"
// #include "Solutions.hpp"

// class fourier
// {
// 	public:
// 	std::vector<std::pair<double,double>> data;
// 	uint size;
// 	double step; 
// };
// /**
//  * @brief 
//  * Class created for passing functions differential equations as parameters, and create states for any given function. 
//  * 
//  * \param n_dim 
//  * Dimension of the differential equation 
//  * \param params
//  * Eventual parameters needed to use ODE.
//  * \param func 
//  * A pointer to the ODE implementation.
//  * \param x0 
//  * Initial conditions of the ODE used in the integration. 
//  * \param solved 
//  * The complete runge-kutta of the ODE
//  * 
//  */
// class function
// {
// 	private:
		
// 		unsigned int n_dim;
// 		double* params;
// 		double* (*func) (std::pair<double*, double>, double*);// F(x,t,param) 
// 	public:
// 		solution* solved = new solution; 
// 		fourier* fourierT = new fourier;
// 		double* operator ()(std::pair<double*, double> x)
// 		{
// 			return func(x, params);
// 		};
// 		// double* operator ()(std::pair<double*, double> x, double* param)
// 		// {
// 		// 	return func(x, param);
// 		// };
// 		function(double* (*fun) (std::pair<double*, double>, double*), std::pair<double*, double> x0, double* param, unsigned int dim);
	
		
// 		// double* runge_kutta_step(std::pair<double*, double> x, double h);
// 		// solution* runge_kutta4th(std::pair<double*, double> x, double tf, double h);
// 		solution* applyRunge_kutta4th(std::pair<double*, double> x0, double t_initial, double t_final, double h);

// 		static double* lorenz_equation(std::pair<double*, double>x, std::vector<double> param);
// 		static double* magenticDipole( std::pair<double*, double> x, double* param);
// 		static double* sistemaAmigos(std::pair<double*, double>x, std::vector<double> param);
// 		std::complex<double> simpsonOneThird(std::complex<double>* f,double t0,double tf, uint size);
// 		std::complex<double> fourierTransform(double frequency);

// 		std::vector<std::pair<std::complex<double>, double>> fourierTransformRange(double inicialFrequency, double finalFrequency, double step, double t0, double tf);
// 		/**
// 		 * @brief Function to find the transform a function in the x domain
// 		 * to the frequency-amplitude domain. The amplitude is the norm squared of the frequencies. 
// 		 * It can take frequencies from f0 to ff, equalilly spaced by step. 
// 		 * 
// 		 * @param step step between frequencies
// 		 * @param f0 initial frequency
// 		 * @param ff final frequency 
// 		 * @return fourier Data type with the transform, size of the transform and step used.
// 		 */
// 		fourier fourierTransformRange(double step, double f0, double ff);
// 		/**
// 		 * @brief Prints the runge kutta already calculated
// 		 * 
// 		 * @param filename name of the file where the runge kutta is to be printed. 
// 		 */
// 		void printRungeKutta(char* filename);
// 		static std::vector<std::pair<std::vector<double>, double>> func_test(unsigned n_points, double t0, double tf);
		
// 		/**
// 		 * @brief function used to find the frequencies in a fourier transform, it works by findind the peaks
// 		 * of the amplitude in the frequency spectrum.
// 		 * 
// 		 * 
// 		 * @param tresholdAmplitude the minimum amplitude of frequencies
// 		 * @return std::vector<std::pair<double,double>> vector of (amplitude, frequency).
// 		 */
// 		std::vector<std::pair<double,double>> findFrquenciesFourier(double tresholdAmplitude);
// 		/**
// 		 * @brief Function to find the maximum value in an array of 2 entries between start and end
// 		 * 
// 		 * @param v vector of 2 entries
// 		 * @param start first entry to be considered
// 		 * @param end last entry to be considered
// 		 * @return std::pair<double,double> pair of maximum values.
// 		 */
//  		static std::pair<double,double> peakSearch(std::vector<std::pair<double,double>> v, int start, int end);
// 		void printFourierTransformRange(char* filename);

// };
