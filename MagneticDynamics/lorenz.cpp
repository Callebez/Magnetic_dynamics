#include"functions.hpp"

class lorenz:function
{
    
public:

    double* param = new double[3]{10.,8./3.,28.};
    function f = function(lorenz::lorenz_equation,3,3,param);

    static double* lorenz_equation(std::pair<double*, double>x, double* param)
    {
        double* x_dot = new double[3];
        x_dot[0] = param[0] * (x.first[1] - x.first[0]);
        x_dot[1] = x.first[0] * (param[2] - x.first[2]) - x.first[1];
        x_dot[2] = x.first[0] * x.first[1] - param[1] * x.first[2];
        
        return x_dot;
    };

};