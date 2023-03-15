#include <iostream>
#include <fstream>
#include<chrono>
#include<cstdlib>
#include "matrices.hpp"
#include "functions.hpp"
#include <ctime>
#include "plotting.h"
#include <algorithm>
#include <cstring>
#include <GL/glut.h>
#include <cmath>
#include <iomanip>
#include <sstream>

#define PI 3.14159265359

// double* Function::test(std::pair<double*, double>x, double* param)
// {
// 	double* y = new double;
// 	*y = 3.0*sin(x.second)+1.5*sin(2.0*x.second)+0.75*sin(3.0*x.second)+0.375*sin(4.0*x.second)+0.01*sin(8.0*x.second);
// 	return y;
// }
std::pair<double,double> largestElement(solution<double>& X)
{
	std::pair<double,double> max = X.data[0];
	for(uint i = 1; i < X.get_n_iterations(); i++)
	{
		if(max.first < X.data[i].first)
		{
			max = X.data[i];
		}
	}
	return max;
}

void simpleFilter(solution<double>& x, std::vector<std::pair<double,double>>& accepted)
{

	for(uint i = 0; i < x.get_n_iterations(); i++)
	{
		if(x.data[i].first/largestElement(x).first >= 0.50)
		{
			accepted.emplace_back(x.data[i]);
		}
	}

}
void ClientCode(std::unique_ptr<Solver>& method, void (*func) (std::pair<double*, double>, double*,double*&), std::pair<double*, double> x0, double* param, double t_initial, double t_final, double h, uint dim)
{
    method->Method(func,x0,param,t_initial, t_final, h, dim);
};

// Window size
const int WIDTH = 640;
const int HEIGHT = 640;

// Structure to store the properties of a circle
struct Circle {
    float x, y, radius,rangle;
};

// float rangle = 0.0f;
// Circle properties
Circle circle1 = { 0.5f, 0.0f, 0.2f,0.0f };
Circle circle2 = { cosf(2.0f*M_PI/3.0f)/2.0f, sinf(2.0f*M_PI/3.0f)/2.0f, 0.2f, 2.0f*M_PI/3.0f};
Circle circle3 = { cosf(4.0f*M_PI/3.0f)/2.0f, sinf(4.0f*M_PI/3.0f)/2.0f, 0.2f, -2.0f*M_PI/3.0f };

// Function to draw a circle
void drawCircle(Circle circle) {
    
	glColor3f(1.0f, 1.0f, 1.0f);
    int segments = 50;
    glBegin(GL_TRIANGLE_FAN);
    glVertex2f(circle.x, circle.y);
    for (int i = 0; i <= segments; i++) {
        float angle = i * M_PI / segments+circle.rangle+M_PI;
        glVertex2f(circle.x + circle.radius * std::cos(angle),
                   circle.y + circle.radius * std::sin(angle));
    }
    glEnd();

    glColor3f(0.5f, 0.5f, 0.5f);
	glBegin(GL_TRIANGLE_FAN);
    	glVertex2f(circle.x, circle.y);
		for (int i = 0; i <= segments; i++) {
			float angle = i * M_PI / segments+ circle.rangle;
			glVertex2f(circle.x + circle.radius * std::cos(angle),
					circle.y + circle.radius * std::sin(angle));
		}
    glEnd();
}


double* param = new double[3]{1.0,1.0,1.0};
TripleMagneticDipole fun = TripleMagneticDipole(param);

int count = 1;
int count_speed;
void updateCircles() {

	count= count + count_speed;
	if(count >fun.solver->rk_int->data.size())
	{
		count = 0.0;
	}

	circle1.rangle = (fun.solver->rk_int->data[count-1].first[0]);
	circle2.rangle = (fun.solver->rk_int->data[count-1].first[1]);
	circle3.rangle = (fun.solver->rk_int->data[count-1].first[2]);
}
void display() {
    glClear(GL_COLOR_BUFFER_BIT);

	glColor3f (1.0, 1.0, 1.0);
	glRasterPos2f(-1.0f, 0.9f); //define position on the screen

	// double pi = 3.14159265359;
	std::stringstream stream;
	stream << std::fixed << std::setprecision(2) << fun.solver->rk_int->data[count].second;
	std::string s = stream.str();

	std::string str = "Tempo de simulação: " + s;
	char *string = new char[str.length()]; 
	strcpy(string, str.c_str()); 
		
	while(*string){
		glutBitmapCharacter(GLUT_BITMAP_8_BY_13, *string++);
	}	

    // Draw circles
    glColor3f(1.0f, 1.0f, 1.0f);
    drawCircle(circle1);
    drawCircle(circle2);
    drawCircle(circle3);
    glutSwapBuffers();
	// delete[] string;
}
// Function to reshape the window
void reshape(int width, int height) {
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0f, 0.0f, 0.0f, 0.0f);

    // gluOrtho2D(-1.0f, 1.0f, -1.0f, 1.0f);
    glMatrixMode(GL_MODELVIEW);
}

// Function called repeatedly to update the scene
void idle() {
	// count=0;
    updateCircles();

	glRasterPos2f(circle1.x, circle1.y);
	// unsigned char string [] = fun.solver->rk_int->data[count-1].second.c_str();
	// std::string str = std::to_string(fun.solver->rk_int->data[count-1].second);
	// const char* string = new char[str.size()];
	// string = str.c_str();
	// int len, i;
	// len = strlen(string);
	// for (i = 0; i < len; i++) 
	// {
	glColor3f(1.0f, 1.0f, 1.0f);
	glutBitmapCharacter(GLUT_BITMAP_8_BY_13	, fun.solver->rk_int->data[count-1].second);
	// }	
    glutPostRedisplay();
	
}
std::string TiraPonto(std::string str)
{
	// std::string ponto= ".";
	// char* ponto = '.';
	// std::string res=str; 
	for(uint i = 0; i <str.length(); i++)
	{
		// if(str[i]==ponto)
		std::cout<<str[i];
	}
	// return str;
}
int main(int argc, char** argv) {


	//////////////////////////////////////

	// // param[2] = 1.0;
	// int time;
	// for(uint i = 0; i < 10; i++)
	// {
	// 	double* param = new double[2];
	// 	param[0] = 0.5; // eps 
	// 	param[1] = 0.5; // omega 
	// 	time = i*30 + 10;
	// 	MagneticDipole fun = MagneticDipole(param);
	// 	std::pair<double*,double > x;
	// 	std::string fileName = "Dipolo_Eps_0_5__Omega_0_5_From_0_"+std::to_string(time);
	// 	std::string fileNameFourier = "Fourier_Dipolo_Eps_0_5__Omega_0_5_From_0_"+std::to_string(time);

	// 	std::string MagDipole = "apresentacao/txts/"+fileName;
	// 	std::string fourier = "apresentacao/txts/"+fileNameFourier;
		
	// 	std::string plotCommand = "plot '" + MagDipole +".txt'"+ " u 1:2 w l lw 3 lc 7 notitle"; 
	// 	std::string plotCommandFourier = "plot '" + fourier +".txt'"+ " w l lw 2 lc 0 notitle"; 
		
	// 	std::string title = "Dipolo magnético na presença de campo magnético \\n constante de frequência de oscilação {/Symbol a}= 0.5";
	// 	std::string titleFourier = "Transformada de fourier para as frequências de oscilação do Dipolo \\n tempo de integração " + std::to_string(time) ;

	// 	double* coord = new double[2];

	// 	coord[0] = 1.0;
	// 	coord[1] = 0.0;
	// 	// coord[2] = -2.094395102393195;
	// 	// coord[3] = 0.00001;
	// 	// coord[4] = 0.0;
	// 	// coord[5] = 0.0;

	// 	x = std::make_pair(coord,0); 
	// 	fun.applySolver<RK4thSolver>(x,0.0,(double) time,1e-3,MagDipole);
		
	// 	plotting::plot2D(fileName,
	// 					plotCommand,
	// 					title,
	// 					"set xlabel \"Tempo característico\" \n set ylabel '{/Symbol q}(t)'");
	// 	fun.applyFourierTransform(-1.5,1.5, 0.01);
	// 	fun.fourier_transform->printSolutionDouble(fourier);
	// 	plotting::plot2D(fileNameFourier,
	// 					plotCommandFourier, 
	// 					titleFourier, 
	// 					"set xlabel \"Frequências\" \n set ylabel 'Amplitudes'");

	// }
	

	double* coord = new double[6];

 	coord[0] = 1.0;
	coord[1] = 0.0;
	coord[2] = -2.094395102393195;
	coord[3] = 2.00001;
	coord[4] = 0.0;
	coord[5] = 0.0;
	std::pair<double*,double > x;
	
	x = std::make_pair(coord,0); 
	fun.applySolver<RK4thSolver>(x,0.0,100.0,1e-3);

	count_speed = 10;

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowSize(640, 480);
	glutCreateWindow("Magnetic Dipoles");

	glutDisplayFunc(display);
	glutIdleFunc(idle);

	glutMainLoop();
	return 0;
}


// int main(void)
// {
// 	auto start = std::chrono::high_resolution_clock::now(); 

// 	std::srand(time(0));
// 	// foo();
// 	// double* x;
// 	// std::unique_ptr<double*> y = foo();
// 	// x =  *(y.get());
// 	// std::cout<<x[0]<<'\n';
// 	// free(x);

// 	// double* x = new double; 
// 	// *x = 5.0;
// 	// double* a = new double; 
// 	// *a = *x; 
// 	// free(x);
// 	// std::cout<<*a<<'\n';

// 	// std::memcpy(x, (*foo().get()), 8);

// 	// std::cout<<(*foo().get())[1];
// 	// std::vector<double> x;
// 	// double* y;
// 	// x = *foo().get();
// 	// y = &(foo().get());
// 	// std::cout<<y;
// 	// std::cout<<foo()->data()[0]<<"\n";
// 	// std::cout<<foo()->data()[1]<<"\n";
// 	// std::cout<<foo()->data()[2]<<"\n";

// 	//////////////////////////////////////
// 	double* param = new double[2];
// 	param[0] = 1.0;
// 	param[1] = 5e-1;
// 	// param[2] = 1.0;

// 	TripleMagneticDipole fun = TripleMagneticDipole(param);
// 	std::pair<double*,double > x;
// 	std::string MagDipole = "tripleDipole.txt";
// 	std::string Tripolefourier = "Tripolefourier.txt";

// 	// for(uint i = 0; i < 10; i++)
// 	// {
// 		double* coord = new double[6];

// 		coord[0] = PI/2.0;
// 		coord[1] = -PI/2.0;
// 		coord[2] = 0.0;
// 		coord[3] = 0.0;
// 		coord[4] = 0.0;
// 		coord[5] = 0.0;
// 		x = std::make_pair(coord,0);
// 		fun.applySolver<RK4thSolver>(x,0.0,500.0,1e-3,MagDipole);
// 		// fun.applyFourierTransform(0,10,0.01);
// 		// fun.fourier_transform->printSolutionDouble(Tripolefourier);
// 		// std::cout<<coord[0]<<'\n';
// 	// }
// 	// std::cout<<((double) rand()) / (RAND_MAX)+1<<"\n";
// 	// coord[2] = 0.0;
// 	// coord[3] = PI;
// 	// coord[4] = 0.00;
// 	// coord[5] = 0.00;

// 	// std::pair<double*,double > x = std::make_pair(coord,0);
// 	// fun.applySolver<RK4thSolver>(x,0.0,1000.0,1e-3,MagDipole);	
// 	// delete [] param;


// 	// free(param);
// 	// delete[] param;
// 	// delete[] coord;
// 	////////////////////////////////////////////////

// 	// fun.applyFourierTransform(0,5,0.01);
// 	// std::string testFourier = "testeFourier.txt";
// 	// fun.fourier_transform->printSolutionDouble(testFourier);
	
// 	// free(param);
// 	// free(coord);

// 	///////////////////////////////////////////////////////////////////////
//  	// std::fstream plot; 
// 	// double step = 0.01;
// 	// double* param = new double[2];
// 	// double* x = new double[2]{1.0,0.0};
// 	// std::vector<std::pair<double,double>> filtered;
// 	// std::pair<double*,double > X = std::make_pair(x,0);
	
// 	// MagneticDipole f = MagneticDipole(param);	
// 	// for(uint i = 0; i < 5; i++)
// 	// {
// 	// 	std::string bif = "bifurcacao" + std::to_string(10+i*20) + ".txt";
// 	// 	plot.open(bif,std::fstream::out);
// 	// 	for(int j = -100; j < 100; j++)
// 	// 		{
// 	// 			param[0] = 1.0;
// 	// 			param[1] = j*step;
// 	// 			f.set_param(param);    
// 	// 			f.applySolver<RK4thSolver>(X,0,10+i*20,1e-3);
				
// 	// 			f.applyFourierTransform(-2,2,0.01);
// 	// 			simpleFilter(*f.fourier_transform,filtered);
			
// 	// 			for(auto k : filtered)
// 	// 			{
// 	// 				plot<<j*step<<" "<<k.first<<" "<<k.second<<"\n";
// 	// 			}
// 	// 			filtered.clear();
// 	// 		}
// 	// 	plot.close();	
// 	// }
// 	// free(x);
// 	// free(param);
// 	////////////////////////////////////////////////
	
// 	// free(solver);;
// 	// double* (*func) (std::pair<double*, double>, double*), std::pair<double*, double> x0, double* param, double t_initial, double t_final, double h, uint dim
// 	// ClienteCode(*solver, Lorenz::lorenz_equation, x, fun.get_param(),0.0,10.0, 0.1, 3);

// 	// ptr = std::unique_ptr<double*>(new double);
	 
	
// 	// x.set_dim(2); 
// 	// for(auto i: *x)
// 	// {
// 	// 	std::cout<<i <<"\n";
// 	// }
// 	// x.time = 5;

// 	// function testFunc = function::test(); 
// 	// double* x = new double[1]; 
// 	// solution<double*>* signal = function::createSignal(testFunc, x,0, 10,0.01,1);
// 	// testFunc.rk_int = signal;
// 	// testFunc.rk_int->t0 = signal->t0;  
// 	// testFunc.rk_int->tf = signal->tf; 
// 	// solution<double>* four = testFunc.applyFourierTransform(0,10,0.001);
// 	// std::string fourierzin = "FourierTest.txt";

// 	// four->printSolutionDouble(fourierzin);



// 	////////////////////////////////////////////////////////////////////////

// 	// double* x = new double[2]{1.0,0.0};
//  	// std::fstream plot; 
// 	// double step = 0.025;
// 	// double* param = new double[2];
// 	// std::vector<std::pair<double,double>> filtered;
// 	// plot.open("bifucacao.txt",std::fstream::out);
	
// 	// std::unique_ptr<function> f = std::make_unique<function>(function::magenticDipole,2,2,param);

// 	// std::pair<double*,double> X = std::make_pair(x,0); 
// 	// for(int j = -100; j < 100; j++)
// 	// {
// 	// 	param[0] = 1.0;
// 	// 	param[1] = j*step;
// 	// 	f->set_param(param);    
// 	// 	f->applyRunge_kutta4th(X,0,100,1e-2);
// 	// 	// solution<double>* four = new solution<double>;
// 	// 	f->applyFourierTransform(-2,2,0.01);
// 	// 	simpleFilter(*f->fourier_transform,filtered);
// 	// 	f->rk_int->data.clear();
// 	// 	f->fourier_transform->data.clear();

// 	// 	plot.open("bifucacao.txt",std::fstream::app);

// 	// 	for(auto k : filtered)
// 	// 	{
// 	// 		plot<<j*step<<" "<<k.first<<" "<<k.second<<"\n";
// 	// 	}
// 	// 	plot.close();

// 	// 	filtered.clear();
// 	// }
// 	// // std::cout<<"salve!\n";
// 	// free(x);
// 	//////////////////////////////////////////////////////////////////


// 	// free(param);
// 	// free(f);

// 	// double* param = new double[2]{1.0, 0.5};
// 	// function* f = new function();
// 	// *f = function(function::magenticDipole,2,2,param);

// 	// solution<double*>* y = new solution<double*>;
// 	// y = f->applyRunge_kutta4th(std::make_pair(x,0),0,50,1e-4);
// 	// solution<double>* four = new solution<double>;
// 	// four = f->applyFourierTransform(-2,2,0.01);
	
// 	// // for(auto i : filtered)
// 	// // {
// 	// // 	std::cout<<i.first<<" "<<i.second<<"\n";
// 	// // }
	
// 	// std::string rk = "Teste.txt";
// 	// std::string fourierzin = "FourierTest.txt";
// 	// y->printSolutionDoublePtr(rk);
// 	// four->printSolutionDouble(fourierzin);
    
// 	// std::fstream plot; 
// 	// plot.open("Dat/Grafico.txt",std::fstream::out);
//     // for(uint i = 0; i < y->n_iterations; i++)
//     // {
//     //     plot<<y->data[i].second <<" "<<y->data[i].first[0]<<" "<<y->data[i].first[1]<<" "<<y->data[i].first[2] <<std::endl; 
//     // }
//     // plot.close();

// 	// std::cout<<d->size<<"\n";
// 	//std::vector<double> x = { 1.0,0.01,0.1};
// 	//std::vector<double> lorenzParams = {10.0,8./3.0,28.0};
// 	//function f = function(function::lorenz_equation, lorenzParams, 3, 3);
// 	//std::vector<std::pair<std::vector<double>, double>> y = f.runge_kutta4th(make_pair(x, 0), 100, 0.001);
// 	//matrix::printMatrixToFile(y, "lorenz.txt");
// 	//plot3D("lorenz", "lorenz", "Sistema de Lorenz", "2:3:4");
	
// 	/////////////////////////////////////////////////////////////
// 	//// Dipolo magn�tico
// 	//std::vector<double> x = { 1.0,0.0 };
// 	//std::vector<double> pendulumParams = { 1.0,2.5};
// 	//
// 	//function magenticDipole = function(function::magenticDipole, pendulumParams, 2);
// 	//
// 	//std::vector<std::pair<std::vector<double>, double>> y = magenticDipole.runge_kutta4th(make_pair(x, 0), 100, 0.001);
// 	//matrix::printMatrixToFile(y, "./Dat/dipole_001_f2.5.txt");
// 	//
// 	//plotting::plot2D("dipole_001_f2.5", "./Dat/dipole_001_f2.5", "Magnetic dipole in an osciling magnectic field with frequency of 2.5", "1:2");
// 	///////////////////////////////////////////////////////////////


// 	// double* x = new double[2];
// 	// x[0] = 1.0;
// 	// x[1] = 0.0;
// 	// double* pendulumParams =new double[2];	
// 	// pendulumParams[0] = 1.0;
// 	// pendulumParams[1] = 0.75;

// 	// function dipole = function(function::magenticDipole, std::make_pair(x,0),pendulumParams,2);
// 	// solution y = *dipole.applyRunge_kutta4th(std::make_pair(x,0),0.0, 100.0*PI,1e-3);
	
// 	// std::fstream plot; 
// 	// plot.open("Dat/Grafico.txt",std::fstream::out);
// 	// for(uint i = 0; i < y.n_iterations; i++)
// 	// {
// 	// 	plot<<y.data[i].second;
// 	// 	for(uint j = 0; j < y.sysDim; j++)
// 	// 	{
// 	// 		plot<<" "<<y.data[i].first[j];
// 	// 	}
// 	// 	plot<<"\n";
// 	// }
// 	// plot.close();
	
// 	// dipole.printRungeKutta("Dat/Grafico.txt");
// 	// dipole.fourierTransformRange(1e-3,0.0,2.0);
// 	// dipole.printFourierTransformRange("Dat/frequenciesTest1.txt");
// 	// std::vector<std::pair<double,double>>frequencies = dipole.findFrquenciesFourier(0.00001);
	
// 	// for(uint i = 0; i<frequencies.size(); i ++)
// 	// {
// 	// 	std::cout<<i << " "<<frequencies[i].first<<" "<<frequencies[i].second<<"\n";
// 	// }

// 	// std::fstream mainFrq;
// 	// mainFrq.open("Dat/mainfrequencies1.txt", std::fstream::out);
	
// 	// std::vector<std::pair<double,double>> F(n);
// 	// F[0].first = std::pow(std::norm(dipole.fourierTransform(0)),2);
// 	// F[1].second = 0.0;

// 	// double step = 1e-3;
// 	// double treshold =0.0001;
// 	// std::vector<std::pair<double,double>> mfreq; 
// 	// std::vector<std::pair<double,double>>v; 

// 	// for(int i = 0 ; i < n; i++)
// 	// {
// 	// 	F[i].first = std::pow(std::norm(dipole.fourierTransform(i*step)),2);
// 	// 	F[i].second = i*step;
// 	// 	frequencies<<F[i].first << " " <<F[i].second <<std::endl;
// 	// }
// 	// frequencies.close();

// 	// for(int i = 0; i < n; i++)
// 	// {
// 	// 	while(F[i].first < treshold && i < n)
// 	// 	{ i++;}
// 	// 	while (F[i].first >= treshold )
// 	// 	{
// 	// 		v.emplace_back(F[i]);
// 	// 		// std::cout<<"i: "<< i <<" " <<v.back().first <<" " <<v.back().second<<"\n";
// 	// 		i++;
// 	// 	}
// 	// 	if(v.size() > 1)
// 	// 	{
// 	// 		mfreq.emplace_back(function::peakSearch(v,0,v.size()));
// 	// 		v.clear();
// 	// 	}
// 	// }

// 	// for(uint i = 0; i < mfreq.size(); i++)
// 	// {
// 	// 	mainFrq<<mfreq[i].first<<" "<< mfreq[i].second<<"\n";
// 	// }
// 	// mainFrq.close();


// 	// double* v = new double[5];
// 	// v[0] = 1.0;
// 	// v[1] = 2.0;
// 	// v[2] = 1.0;
// 	// v[3] = 3.0;
// 	// v[4] = 1.0;
// 	// std::cout<<"Pico "<<peakSearch(v,0,5);
// 	// std::fstream output;
// 	// output.open("Dat/teste1.txt", std::fstream::out);
// 	// for(uint i = 0; i < dipole.solved.size; i++)
// 	// {
// 	// 	output<<rk.data[i].first[0]<<", "<<rk.data[i].first[1]<<", "<<rk.data[i].second<<std::endl;
// 	// }
// 	// output.close();


// 	auto stop = std::chrono::high_resolution_clock::now();
// 	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
	
// 	std::cout << "\n Time taken by the program: " << duration.count() << " miliseconds" << "\n";

// 	return 0;
// }