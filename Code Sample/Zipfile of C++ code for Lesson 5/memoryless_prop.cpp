// An illusration of the memoryless property of exponential distributions
// Written by Prof. Sreenivas for IE523: Financial Computing
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <random>
#include <chrono>

using namespace std;

unsigned seed = (unsigned) std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator;

double get_exp(double lambda)
{
    // http://www.cplusplus.com/reference/random/exponential_distribution/
    std::exponential_distribution<double> distribution(lambda);
    double number = distribution(generator);
    return (number);
}

int main (int argc, char* argv[])
{
	double lambda, cut_off_value, x;
	int no_of_trials, no_greater_than_cut_off, count[100];
	
	sscanf (argv[1], "%lf", &lambda);
	sscanf (argv[2], "%d", &no_of_trials);
	sscanf (argv[3], "%lf", &cut_off_value);
	ofstream outf(argv[4]);
	
	no_greater_than_cut_off = 0;
	for (int i = 0; i < 100; i++)
		count[i] = 0;
	
	for (int i = 0; i < no_of_trials; i++) {
		x = get_exp(lambda);
		if (x > cut_off_value) {
			no_greater_than_cut_off++;
			for (int j = 0; j < 100; j++)
				if ( (x - cut_off_value) < ((float) j/10))
					count[j]++;
		}
	}
	
    cout << "Testing the Memoryless Property of Exponential Distributions" << endl;
    cout << "Rate (i.e. Lambda)                       = " << lambda << endl;
    cout << "Number of samples >= cut-off-value of " << cut_off_value << " = " << no_greater_than_cut_off << endl;
    if (no_greater_than_cut_off > 0)
    {
        for (int j = 0; j < 100; j++)
            outf << ((double) j/10) << ", " <<
            ((double) count[j])/((double) no_greater_than_cut_off) << ", " <<
            (1.0-exp(-1.0*lambda*((double) j/10))) << endl;
    }
    else
    {
        cout << "Insufficient #samples >= " << cut_off_value << " to form an empirical-distribution " << endl;
        cout << "Reduce the cut-off value and try again!" << endl << endl;
    }

}


