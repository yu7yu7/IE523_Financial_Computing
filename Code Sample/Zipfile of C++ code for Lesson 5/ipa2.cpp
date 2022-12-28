// Illustration of Infinitesimal Perturbation Analysis via an
// Example.  Simulation Program #2.
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
	float x, y, theta, derivative_y_wrt_theta; 
	sscanf (argv[1], "%d", &no_of_trials); 
	sscanf (argv[2], "%f", &theta); 
	cout << "\nIllustration of IPA..."; 
	cout << "\nTheta = " << theta << endl; 
	
	y = 0.0; 
	y = derivative_y_wrt_theta = 0.0; 
	for (int i = 0; i < no_of_trials; i++) { 
		x = get_uniform(); 
		y = y + (theta * x); 
		derivative_y_wrt_theta = derivative_y_wrt_theta + x; 
	} 
	cout << "Average value of y for " << no_of_trials << " trials = " 
	<<	y/((float) no_of_trials) << endl;
	cout << "\nAverage derivative of y wrt theta for " << no_of_trials << 
	" trials = " << derivative_y_wrt_theta/((float) no_of_trials) << endl; 
}

	
