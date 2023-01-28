#pragma once 
#include "RungeKutta.hpp"
#include "Solutions.hpp"
#include "fourier.hpp"


class Function
{
private:
    uint x_dim = 0;
    uint f_dim = 0;
    double* (*func) (std::pair<double*, double>, double*) = NULL;
protected:
    // std::unique_ptr<Solver> solver; 
public:
    double* params = NULL;
	
    inline uint get_x_dim(){return x_dim;};
	inline uint get_f_dim(){return f_dim;};
    inline double* get_param(){return params;};
    // double* (*func) (std::pair<double*, double>, double*) get_func(){return func;};
	void set_x_dim(uint new_x_dim){x_dim =  new_x_dim;};
	void set_f_dim(uint new_f_dim){f_dim =  new_f_dim;};
    void set_func(double* (*fun) (std::pair<double*, double>, double*)){func = fun;};
    void set_param(double* param){params = param;};
    // void set_param(double value, uint pos){params[pos] = value;};
    Function(){;};
    Function(double* (*func_) (std::pair<double*, double>, double*),uint x_dim_,uint f_dim_, double* param_)
    {
        x_dim = x_dim_;
        f_dim = f_dim_;
        params = param_;
        func = func_;
    };

    double* operator()(std::pair<double*, double> x)
    {
        return func(x,params);
    };
    void print()
    {
        std::cout<<"params: ";
        for(uint i = 0; i < 3; i++)
        {
            std::cout<<params[i]<<" ";
        }
        std::cout<<std::endl; 
        std::cout<<"dim x: "<< get_x_dim()<<", dim f: "<<get_f_dim()<<std::endl;
        std::cout<<"func: "<<&func;
    };
    
    // void applySolver(Solver& solverMethod, double t0, double tf, double step_inicial);
    // virtual void runge_Kutta(std::pair<double*, double> x0, double t_initial, double t_final, double h);
    // virtual void Fourier_Transform(double initialFrequency, double finalFrequency, double frquencyStep);
    
};
class Lorenz : public Function
{
public: 
    /** Constructor **/
    Lorenz(double* param)
    {
        // std::cout<<"salve\n";
        set_param(param);
        set_func(lorenz_equation);
        set_x_dim(3);
        set_f_dim(3);
    };
    
    /**Differential equation for the Lorenz system **/    
    static double* lorenz_equation(std::pair<double*, double>x, double* param);

};
// /**
//  * @brief Base function for the creation of concret functions. 
//  * 
//  * @param x_dim Dimension of the domain;
//  * @param f_dim Dimension of the image;
//  * @param params Parameters of the function;
//  * @param func pointer to the method where the function is implemented;
//  * @param rk_int pointer to the solution of the system of ODEs;
//  * @param fourier_transform pointer to the solution of the fourier transform; 
//  * 
//  */
// class function
// {
// private:
//     uint x_dim = 0; 
//     uint f_dim = 0;
//     double* params = NULL; 
// 	double* (*func) (std::pair<double*, double>, double*) = NULL;// F(x,t,param) 
// public:

//     std::unique_ptr<solution<double*>> rk_int;
//     std::unique_ptr<solution<double>> fourier_transform;
//     static double* test(std::pair<double*, double>x, double* param);
//     function(){;};
//     static inline function test()
//     {
//         double* param = new double;
//         function testf = function(test,1,1,param);
//         return testf;
//     };
//     static inline function lorenz()
//     {
//         double* param = new double[3]{10.,8./3.,28.};
//         function lorenz = function(lorenz_equation, 3,3,param);
//         return lorenz;
//     };
//     static inline function magneticDipole()
//     {
//         double* param = new double[2]{1.,2.5};
//         function magneticDipole = function(magenticDipole, 2,2,param);
//         return magneticDipole;
//     };

//    /**
//     * @brief Construct a new function object
//     * 
//     * @param fun Method where the function is implemented;
//     * @param xdim Dimension of the domain;
//     * @param fdim Dimension of the image;
//     * @param param  Parameters of the function;
//     */
//     function(double* (*fun) (std::pair<double*, double>, double*),uint xdim,uint fdim, double* param)
//     {
//         x_dim = xdim;
//         f_dim = fdim;
//         params = param;
//         func = fun;
//     };
//     // ~function()
//     // {
//     //     std::cout<<"The destructor was called!\n";
//     //     delete rk_int;
//     //     delete fourier_transform;
//     // }
//     /**
//      * @brief Overloads the operator () so that the pointer func can be used to evaluete the function.
//      * 
//      * @param x Point where the function is being evalueted.
//      * @return double* f(x).
//      */
//     double* operator()(std::pair<double*, double> x)
//     {
//         return func(x,params);
//     };

// 	inline uint get_x_dim(){return x_dim;};
// 	inline uint get_f_dim(){return f_dim;};
//     inline double* get_param(){return params;};

// 	inline void set_x_dim(uint new_x_dim){x_dim =  new_x_dim;};
// 	inline void set_f_dim(uint new_f_dim){f_dim =  new_f_dim;};
//     inline void set_func(double* (*fun) (std::pair<double*, double>, double*)){func = fun;};
//     inline void set_param(double* param){params = param;};
//     inline void set_param(double value, uint pos){params[pos] = value;};
// /**
//  * @brief Create a Signal object
//  * 
//  * @param signalFunc Function used for creating the signal 
//  * @param params_ parameters of the signal function
//  * @param t0_ initial time of the signal
//  * @param tf_ final time of the signal 
//  * @param step_ step of signal sampling 
//  * @param sysDim_ dimension of the signal
//  * @return solution<double*>* signal
//  */
//     static solution<double*>*  createSignal(function signalFunc, double* params_, double t0_, double tf_, double step_, uint sysDim_);

//     /**
//      * @brief Function to calculate the Runge Kutta 4th order of a given system of EDOs
//      * 
//      * @param x0 Iinitial position;
//      * @param t_initial Initial time of integration;
//      * @param t_final Final time of integration;
//      * @param h Step used in the integration; 
//      * @return solution* List of coordinates (X, t) for the integration.
//      */
//     void applyRunge_kutta4th(std::pair<double*, double> x0, double t_initial, double t_final, double h);
//     /**
//      * @brief Function to apply the Fourier Transform of a function. It only works in function objects! 
//      * 
//      * @param initialFrequency Initial frequency in the Fourier transform, lower limit of the integral;
//      * @param finalFrequency Final frequency in the Fourier transform, upper limit of the integral;
//      * @param frquencyStep Step of the frequency. Step of the integration; 
//      * @return solution<double>* Fourier transform of a function from the initial frequency to the final frequency.
//      */
//     void applyFourierTransform(double initialFrequency, double finalFrequency, double frquencyStep);
//     /**
//      * @brief Implementation of the lorenz system (https://en.wikipedia.org/wiki/Lorenz_system)
//      * 
//      * @param x Coordinates;
//      * @param param Paramaters;
//      * @return double* F(x).
//      */
//     static double* lorenz_equation(std::pair<double*, double>x, double* param);
//     static double* magenticDipole( std::pair<double*, double> x, double* param);
     
// };

