#include "plotting.h"

// void plotting::plot3D(const char* outputFile, const char* inputFile, const char* graphTitle, const char* usingLines)
// {
// 	FILE* gnupipe = popen("gnuplot -persist", "w");
// 	if (gnupipe)
// 	{
// 		fprintf(gnupipe, "set term pngcairo enhanced size 1080,720 font \"Times-New-Roman,16\"\n\n");
// 		fprintf(gnupipe, "set title \"%s\"\n\n ", graphTitle);
// 		//fprintf_s(gnupipe, "set title font \"Helvetica, 16\"\n\n");
// 		fprintf(gnupipe, "set tics font \'Times-New-Roman,16\' \n\n");
// 		fprintf(gnupipe, "set border lw 3 \n \n");
// 		fprintf(gnupipe, "set key opaque \n\n");
// 		fprintf(gnupipe, "set output \'%s.png\'\n", outputFile);
// 		fprintf(gnupipe, "splot \'%s.txt\' u 2:3:4 w lines lw 0.5 lc 7 title \'%s\'\n", plotCommand, usingLines);
// 	}
// }

void  plotting::plot2D(std::string outputFile, std::string plotCommand, std::string graphTitle, const char* optionalCommands)
{
	FILE* gnupipe = popen("gnuplot -persist", "w");

	if (gnupipe)
	{
		fprintf(gnupipe, "set encoding utf8 \n\n set terminal png enhanced size 1080,720 font \"Times-New-Roman,16\"\n\n");
		fprintf(gnupipe, "set title font \'Times-New-Roman,24\' \n\n");
		fprintf(gnupipe, "set title \"%s\" \n\n ", graphTitle.c_str());
		
		fprintf(gnupipe, "set tics font \'Times-New-Roman,20\' \n\n");
		fprintf(gnupipe, "set xlabel font \'Times-New-Roman,20\' \n\n");
		fprintf(gnupipe, "set ylabel font \'Times-New-Roman,20\' \n\n");
		fprintf(gnupipe, "set border lw 2 \n \n");
		fprintf(gnupipe, "set grid \n\n");
		// fprintf(gnupipe, "set key noautotitle \n\n");
		fprintf(gnupipe, "set output \'./Omega/images/%s.png\' \n\n",outputFile.c_str());
		fprintf(gnupipe, "%s\n", optionalCommands); 
		fprintf(gnupipe, "%s\n", plotCommand.c_str());
	}

}
