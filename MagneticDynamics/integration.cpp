#include "integration.hpp"
std::complex<double> integration::simpsonOneThird(std::complex<double>* f,double t0,double tf, uint size)
{
    unsigned int n;
    if (size % 2 != 0) { n = size - 1; }
    else { n = size;}   
    std::complex<double> h = (tf - t0) / (double)n;
    std::complex<double> sum = {0.0, 0.0};
    sum += f[0];
    sum += f[n];
    for (unsigned int i = 1; i < uint((n-1)/2)+1; i++)
    {
        sum += f[2 * i - 1] * 4.0 ;
        sum += f[2 * i]* 2.0;
    }
    return sum*h/3.0;
}