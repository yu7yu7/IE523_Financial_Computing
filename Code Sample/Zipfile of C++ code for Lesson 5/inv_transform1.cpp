// Illustration of the Inverse Transform Method
// IID Exponential RV generation
// In practice, you will just use the C++ STL directly
// for Exponential RVs.  This is just an illustration

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <random>
#include <chrono>

using namespace std;

// cf http://www.cplusplus.com/reference/random/uniform_real_distribution/operator()/
// If you want to set a seed -- do it only after debug phase is completed
// otherwise errors will not be repeatable.
//unsigned seed = (unsigned) std::chrono::system_clock::now().time_since_epoch().count();
//default_random_engine generator (seed);

default_random_engine generator;

double get_uniform()
{
    std::uniform_real_distribution <double> distribution(0.0,1.0);
    double number = distribution(generator);
    return (number);
}

double get_exp(double lambda)
{
	return ((-1.0/lambda)*log(1.0-get_uniform()));
}

int main (int argc, char* argv[])
{
	double lambda;
	int no_of_trials, count[100];
	
	sscanf (argv[1], "%lf", &lambda);
	sscanf (argv[2], "%d", &no_of_trials);
	ofstream outf(argv[3]);
    
	for (int i = 0; i < 100; i++)
		count[i] = 0;
	
	for (int i = 0; i < no_of_trials; i++) {
		double y = get_exp(lambda);
		for (int j = 0; j < 100; j++)
			if (y < ((double) j/10))
				count[j]++;
	}
	
	for (int j = 0; j < 100; j++)
		outf << ((double) j/10) <<
        ", " << ((double) count[j])/((double) no_of_trials) << ", " <<
        (1.0-exp(-1.0*lambda*((double) j/10))) << endl;

}


