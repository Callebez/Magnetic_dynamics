// #include "fourier.hpp"
// // void fourier::fourierFrequencySpectrumAbsoluteValue(solution<double *> *signal, double initialFrequency, double finalFrenquency, double frequencyStep)

// void fourier::fourierFrequencySpectrumAbsoluteValue(std::unique_ptr<solution<double*>>& signal, double initialFrequency, double finalFrequency, double frequencyStep)
// {
//     fourierTransform = std::make_unique<solution<double>>();
//     fourierTransform->step = frequencyStep;
//     fourierTransform->t0 = initialFrequency;
//     fourierTransform->tf = finalFrequency;
//     fourierTransform->sysDim = 1;  
//     fourierTransform->n_iterations = (uint) ((finalFrequency-initialFrequency)/frequencyStep);

//     fourierTransform->data = std::vector<std::pair<double,double>>(fourierTransform->n_iterations);
 
//     for(uint i = 0; i < fourierTransform->n_iterations; i++)
//     {      
//         fourierTransform->data[i].first = std::pow(std::abs(fourierFreqency(signal, initialFrequency + i*frequencyStep)),2);
//         fourierTransform->data[i].second = initialFrequency + i*frequencyStep;
//     }

// }
// std::complex<double> fourier::fourierFreqency(std::unique_ptr<solution<double*>>& signal, double frequency)
// {
//     std::vector<std::complex<double>> fexp (signal->n_iterations);
//     // std::complex<double>* fexp = new std::complex<double> [signal->n_iterations];
//     std::complex<double> transform ;
//     for (unsigned int i = 0; i < signal->n_iterations; i  ++)
//     {
//         fexp[i] = std::exp(std::complex<double>(0,-frequency*signal->data[i].second))* signal->data[i].first[0];

//     }
//     // std::cout<<signal->tf<<"\n";

//     transform = 2.0/(signal->tf-signal->t0)*integration::simpsonOneThird(fexp, signal->t0, signal->tf ,signal->n_iterations);
//     // std::cout<<transform<<"\n";

//     return transform;
// }
// // void fourier<std::complex<double>>::fourierFrequencySpectrum(solution<double*>* signal, double initialFrequency, double finalFrequency, double frequencyStep)
// // {
// //     std::cout<<"estou aqui bb\n";
// //     fourierTransform->step = frequencyStep;
// //     fourierTransform->t0 = initialFrequency;
// //     fourierTransform->tf = finalFrequency;
// //     fourierTransform->sysDim = 1; 
// //     fourierTransform->n_iterations = (uint) ((finalFrequency-initialFrequency)/frequencyStep);
// //     for(uint i = 0; i < signal->n_iterations; i++)
// //     {      
// //         fourierTransform->data[i].first[0] = fourierFreqency(signal, initialFrequency + i*frequencyStep);
// //         fourierTransform->data[i].second = initialFrequency + i*frequencyStep;
// //     }
// // }

