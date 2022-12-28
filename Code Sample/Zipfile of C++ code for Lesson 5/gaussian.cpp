// Unit Normal variate Generator using the Box-Muller Transform
// Written by Prof. Sreenivas for IE523: Financial Computing
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <random>
#include <chrono>
using namespace std;

#define GAUSSIAN_DENSITY(x) 0.398942280375387*exp(-x*x/2.0) 

unsigned seed = (unsigned) std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator(seed);

double get_uniform()
{
    // http://www.cplusplus.com/reference/random/exponential_distribution/
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    double number = distribution(generator);
    return (number);
}

double get_gaussian()
{
    return (sqrt(-2.0*log(get_uniform()))*cos(6.283185307999998*get_uniform()));
}

int main (int argc, char* argv[])
{
	float y, cdf_gaussian[100];
	int no_of_trials, count[100];
	
	sscanf (argv[1], "%d", &no_of_trials);
	ofstream cdf_comparison_file(argv[2]);
	ofstream pdf_comparison_file(argv[3]);
	
	for (int i = 0; i < 100; i++) {
		count[i] = 0;
		cdf_gaussian[i] = 0.0;
	}
	
	for (int i = 0; i < no_of_trials; i++) {
		y = get_gaussian();
		for (int j = 0; j < 100; j++)
			if ( (y >= ((float) (j-51)/10)) && (y < ((float) (j-50)/10)) )
				count[j]++;
	}
	
	// CDF of gaussian
	for (int j = 1; j < 100; j++)
		cdf_gaussian[j] = cdf_gaussian[j-1] + 
		(0.1*GAUSSIAN_DENSITY((float) (j-50)/10));
	
	int sum = 0;
	for (int j = 0; j < 100; j++) {
		sum += count[j];
		cdf_comparison_file << ((float) (j-50)/10) << ", " << ((float) sum/no_of_trials) << ", " 
		<< cdf_gaussian[j] << endl;
		pdf_comparison_file << ((float) (j-50)/10) << ", " << 
		((float) count[j]/no_of_trials) << ", " << 
		(0.1*GAUSSIAN_DENSITY((float) (j-50)/10)) << endl;
	}
}
