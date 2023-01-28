// #pragma once
// #include "Solutions.hpp"
// #include "integration.hpp"


// class fourier
// {
// // private:    
// public: 
//     std::unique_ptr<solution<double>> fourierTransform;
//     /**
//      * @brief Integral of a function time a complex exponential. 
//      * 
//      * This function is only used in Fourier Frequency Spectrum. 
//      * 
//      * @param signal value of the function;
//      * @param frequency exponent of the complex exponential 
//      * @return std::complex<double> value of the integral
//      */
//     std::complex<double> fourierFreqency(std::unique_ptr<solution<double*>>& signal, double frequency);

//     /**
//      * @brief Fourier Transform of a signal of function. 
//      * 
//      * @param signal Signal or array of function evaluated in a domain;
//      * @param initialFrequency Lower limit of the Fourier Transform;
//      * @param finalFrenquency Upper limite of the Fourier Transform;
//      * @param frequencyStep Step of the integration;
//      */
//     void fourierFrequencySpectrumAbsoluteValue(std::unique_ptr<solution<double*>>& signal, double initialFrequency, double finalFrequency, double frequencyStep);
//     std::unique_ptr<solution<double>> get_FrequencySpectrum(){return std::move(fourierTransform);};
//     // ~fourier()
//     // {
//     //     fourierTransform->data.clear();
//     //     // delete fourierTransform->data; 
//     //     // delete fourierTransform;
//     //     // delete this;
//     //     // std::cout<<"the fourier transform destructor was called!\n";
//     // }
// };