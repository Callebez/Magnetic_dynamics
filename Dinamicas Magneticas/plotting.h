#pragma once
#include<cstdlib>
#include<fstream>
class plotting
{
public:
	static void plot3D(const char* outputFile, const char* inputFile, const char* graphTitle, const char* usingLines);
	static void plot2D(const char* outputFile, const char* inputFile, const char* graphTitle, const char* usingLines);
};

