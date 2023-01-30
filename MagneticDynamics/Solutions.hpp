#pragma once
#include"matrices.hpp"
#include <memory>


template <class T>
class solution
{
public:
	std::vector<std::pair<T,double>> data;
	
	double t0 ;
    double tf ; 
    double step ; 
    double* params;
    uint n_iterations; 
    uint sysDim; 

	void printSolutionDoublePtr(std::string& filename)
	{
	    std::fstream plot; 
		plot.open(filename,std::fstream::out);
	
		for(uint i = 0; i < data.size()-1; i++)
		{
			plot<<data[i].second;

			for(uint j = 0; j < sysDim;j++)
			{
				plot<<" "<<data[i].first[j];
			}
			plot<<"\n";
		}

		plot.close();
	};
	void printSolutionDouble(std::string& filename)
	{
	    std::fstream plot; 
		plot.open(filename,std::fstream::out);
	
		for(uint i = 0; i < n_iterations; i++)
		{
			plot<<data[i].second;
			plot<<" "<<data[i].first;
			plot<<"\n";
		}
		plot.close();
	};
	inline double get_t0(){return t0;};
	inline double get_tf(){return tf;};
	inline double get_step(){return step;};
	inline double* get_params(){return params;};
	inline uint get_n_iterations(){return n_iterations;};
	inline uint get_sysDim(){return sysDim;};

	inline void set_t0(double new_t0){t0 = new_t0;};
	inline void set_tf(double new_tf){tf = new_tf;};
	inline void set_step(double new_step){ step = new_step;};
	inline void set_params(double* new_param){params = new_param;};
	inline void set_n_iterations(uint new_n_iterations ){ n_iterations = new_n_iterations;};
	inline void set_sysDim(uint new_sysDim){sysDim = new_sysDim;};
};