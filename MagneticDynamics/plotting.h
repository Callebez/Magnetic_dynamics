#pragma once
#include<cstdlib>
#include<fstream>
#include<iostream>
class plotting
{
public:
	static void plot3D(const char* outputFile, const char* inputFile, const char* graphTitle, const char* usingLines);
	static void plot2D(std::string outputFile, std::string inputFile, std::string graphTitle, const char* optionalCommands);

};

