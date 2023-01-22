#pragma once 
#include "RungeKutta.hpp"
#include "Solutions.hpp"
#include "fourier.hpp"

/**
 * @brief Base function for the creation of concret functions. 
 * 
 * @param x_dim Dimension of the domain;
 * @param f_dim Dimension of the image;
 * @param params Parameters of the function;
 * @param func pointer to the method where the function is implemented;
 * @param rk_int pointer to the solution of the system of ODEs;
 * @param fourier_transform pointer to the solution of the fourier transform; 
 * 
 */
class function
{
private:
    uint x_dim = 0; 
    uint f_dim = 0;
    double* params = NULL; 
	double* (*func) (std::pair<double*, double>, double*) = NULL;// F(x,t,param) 
public:

    solution<double*>* rk_int = new solution<double*>;
    solution<double>* fourier_transform = new solution<double>;
    static double* test(std::pair<double*, double>x, double* param);
    function(){;};
    static inline function test()
    {
        double* param = new double;
        function testf = function(test,1,1,param);
        return testf;
    };
    static inline function lorenz()
    {
        double* param = new double[3]{10.,8./3.,28.};
        function lorenz = function(lorenz_equation, 3,3,param);
        return lorenz;
    };
    static inline function magneticDipole()
    {
        double* param = new double[2]{1.,2.5};
        function magneticDipole = function(magenticDipole, 2,2,param);
        return magneticDipole;
    };

   /**
    * @brief Construct a new function object
    * 
    * @param fun Method where the function is implemented;
    * @param xdim Dimension of the domain;
    * @param fdim Dimension of the image;
    * @param param  Parameters of the function;
    */
    function(double* (*fun) (std::pair<double*, double>, double*),uint xdim,uint fdim, double* param)
    {
        x_dim = xdim;
        f_dim = fdim;
        params = param;
        func = fun;
    };
    // ~function()
    // {
    //     delete rk_int;
    //     delete fourier_transform;
    // }
    /**
     * @brief Overloads the operator () so that the pointer func can be used to evaluete the function.
     * 
     * @param x Point where the function is being evalueted.
     * @return double* f(x).
     */
    double* operator()(std::pair<double*, double> x)
    {
        return func(x,params);
    };

	inline uint get_x_dim(){return x_dim;};
	inline uint get_f_dim(){return f_dim;};
    inline double* get_param(){return params;};

	inline void set_x_dim(uint new_x_dim){x_dim =  new_x_dim;};
	inline void set_f_dim(uint new_f_dim){f_dim =  new_f_dim;};
    inline void set_func(double* (*fun) (std::pair<double*, double>, double*)){func = fun;};
    inline void set_param(double* param){params = param;};
    inline void set_param(double value, uint pos){params[pos] = value;};
/**
 * @brief Create a Signal object
 * 
 * @param signalFunc Function used for creating the signal 
 * @param params_ parameters of the signal function
 * @param t0_ initial time of the signal
 * @param tf_ final time of the signal 
 * @param step_ step of signal sampling 
 * @param sysDim_ dimension of the signal
 * @return solution<double*>* signal
 */
    static solution<double*>*  createSignal(function signalFunc, double* params_, double t0_, double tf_, double step_, uint sysDim_);

    /**
     * @brief Function to calculate the Runge Kutta 4th order of a given system of EDOs
     * 
     * @param x0 Iinitial position;
     * @param t_initial Initial time of integration;
     * @param t_final Final time of integration;
     * @param h Step used in the integration; 
     * @return solution* List of coordinates (X, t) for the integration.
     */
    solution<double*>* applyRunge_kutta4th(std::pair<double*, double> x0, double t_initial, double t_final, double h);
    /**
     * @brief Function to apply the Fourier Transform of a function. It only works in function objects! 
     * 
     * @param initialFrequency Initial frequency in the Fourier transform, lower limit of the integral;
     * @param finalFrequency Final frequency in the Fourier transform, upper limit of the integral;
     * @param frquencyStep Step of the frequency. Step of the integration; 
     * @return solution<double>* Fourier transform of a function from the initial frequency to the final frequency.
     */
    solution<double>* applyFourierTransform(double initialFrequency, double finalFrequency, double frquencyStep);
    /**
     * @brief Implementation of the lorenz system (https://en.wikipedia.org/wiki/Lorenz_system)
     * 
     * @param x Coordinates;
     * @param param Paramaters;
     * @return double* F(x).
     */
    static double* lorenz_equation(std::pair<double*, double>x, double* param);
    static double* magenticDipole( std::pair<double*, double> x, double* param);
     
};

