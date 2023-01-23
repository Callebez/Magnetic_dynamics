#pragma once
#include<vector>
#include "matrices.hpp"
#include "functions.hpp"
#include "Solutions.hpp"


/**
 * @brief Class of runge kutta methods. 
 * This class is only used inside the concret class of function, 
 * it solves differential equation given as arguments in concret class.
 * 
 * After a runge kutta method is applied, the solution is stored in "solution"
 * pointer, as a list consisting of the coordinates (x_1,x_2,...,x_n,t) = (X,t). 
 * 
 * \param t0 Initial time of integration;
 * \param tf Final time of integration;
 * \param step Initial step used for the integrations;
 * \param n_iterations size of the integration process. Number of times that a runge kutta method was applied;
 * \param dim Number of differential equations in the system;
 * \param solution List of coordinates (X,t) with the same size of n_iterations;
 * 
 */
class rungeKutta
{
private:
public:
    solution<double*> rk_int;  

    /**
     * @brief Function to perform the a step of runge kutta integration, 
     * \f$ y_{n+1} = y_n + \frac{1}{6}(k_1+k_2+k_3+k_4)*h \f$
     * where \f$k_1 = f(y_n,t_n)\f$, \f$k_2 = f(k_1 + h/2,t_n + h/2)\f$,
     * \f$k_3 = f(k_2 + h/2,t_n + h/2)\f$ \f$k_4 = f(k_3 + h,t_n + h)\f$.
     * 
     * @param func System of differential equations to be solved;
     * @param x Point in which the EDO will be solved;
     * @return double* Coordinates of the next point. 
     */
    double* runge_kutta_step(double* (*func) (std::pair<double*, double>, double*), std::pair<double*, double> x);
    
    /**
     * @brief Runge-Kutta 4 order applied to a function in a given range of time. 
     * Uses the "runge_kuuta_step" function to solve the ODE from t0 to tf. 
     * Only meant to be used in the concret function Class.
     * 
     * @param func System of ODEs F(X,t); 
     * @param x0 Initial point of integration;
     * @param param Parameters of the system of equations;
     * @param t_initial Initial time of integration;
     * @param t_final Final time of integration;
     * @param h Step;
     * @param dim Dimension of the system of ODEs;
     */
    void runge_kutta4th(double* (*func) (std::pair<double*, double>, double*), std::pair<double*, double> x0, double* param, double t_initial, double t_final, double h, uint dim);
    solution<double*> getSolution(){return *rk_int;};
    

};