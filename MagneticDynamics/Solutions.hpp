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
		std::string file = filename+".dat";
		plot.open(file,std::fstream::out);
	
		for(uint i = 0; i < data.size(); i++)
		{
			plot<<data[i].second;

			for(uint j = 0; j < sysDim;j++)
			{
				plot<<" "<<data[i].first[j];
			}
			plot<<"\n";
		}
		plot<<"\n\n\n";

		plot.close();
	};
	void printSolutionDouble(std::string& filename)
	{
	    std::fstream plot; 
		std::string file = filename+".txt";

		plot.open(file,std::fstream::out);
	
		for(uint i = 0; i < n_iterations; i++)
		{
			plot<<data[i].second;
			plot<<" "<<data[i].first;
			plot<<"\n";
		}
		plot.close();

	};
	void filter(std::vector<std::pair<double,double>>& res, double frequencyTreshold)
	{
		std::vector<std::pair<double,double>> aux; 

		for(uint i = 1; i < data.size()-1; i++)
		{
			if((data[i-1].first  < data[i].first )&& (data[i].first  > data[i+1].first ))
			{
				aux.emplace_back(data[i]);
			}
		}
		std::pair<double,double> max = aux[0];
		for(uint i = 1; i < aux.size(); i++)
		{
			if(max.first < aux[i].first)
			{
				max = aux[i];
			}
		}
		for(uint i = 0; i < aux.size(); i++)
		{
			if(aux[i].first/max.first > frequencyTreshold)
			{
				res.emplace_back(aux[i]);
			}
		}

	};
	~solution()=default; 
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