#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <math.h>
#include <algorithm>
#include <msclr\marshal_cppstd.h> //for converting System::String^ to string 
#include <windows.h> //change file directory
#include <tchar.h>
#include <iomanip>
#include <complex>
#include <limits>
#include <mutex>
#include <thread>
#include <omp.h>
#include <functional>
#include <process.h> //run external process
//Sth to external executable acces
#include <shlobj.h>
#include <shlwapi.h>
#include <objbase.h>

class Hamiltonian {
private:
	double omega;
	double QW_position;
	double alfa;

	double x_1; 
	double x_2; 
	double dx;
	int N; // number of points on axis
	
	double E_min, E_max;

public:
	std::vector<std::vector<double>> H;
	std::vector<double> eigenVal;

	Hamiltonian(double omega, double QW_pos, double x_1, double x_2, double dx);
	~Hamiltonian();

	void Build_Hamiltonian(std::vector<double> potential);
	void Strum_Ehrlich_Diagonalization(int num_of_eigenVal);
	void Sturm_Diagonalization(int num_of_eigenVal);

	double Sturm_sequence(double x, int k);
	double Sturm_sequence_derivative(double x, int k);

	std::vector<double> Sturm_sequence_sign_change(double x);
	std::vector<double> find_set_containing_eigenVal(int k);

	std::vector<double> Calculate_eigenvector_INVERSE_POWER(int n); //calculates the n-th eigenvector
	std::vector<double> Calculate_eigenvector(int n); //calculates the n-th eigenvector
	std::vector<double> Calculate_eigenvector_Numerov(std::vector<double> potential, int n); //calculates the n-th eigenvector with NUmerov method
	std::vector<double> Analitycal_eigenvector(int n); //calculates thwe n-th analitycal eigenvector
};

long long int factorial(int n);