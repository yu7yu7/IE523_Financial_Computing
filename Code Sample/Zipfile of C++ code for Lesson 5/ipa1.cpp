// Illustration of Infinitesimal Perturbation Analysis via an
// Example.  Simulation Program #1.
// Written by Prof. Sreenivas for IE523: Financial Computing
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
using namespace std;

float get_uniform()
{
	return (((float) random())/(pow(2.0, 31.0)-1.0));
}

int main (int argc, char* argv[])
{
	int no_of_trials; 
	float x, y, theta; 
	sscanf (argv[1], "%d", &no_of_trials); 
	sscanf (argv[2], "%f", &theta); 
	y = 0.0; 
	for (int i = 0; i < no_of_trials; i++) { 
		x = get_uniform(); 
		y = y + (theta * x); 
	} 
	cout << "\nAverage value of y for " << no_of_trials << " trials = " 
	<<	y/((float) no_of_trials) << endl;
}

	
