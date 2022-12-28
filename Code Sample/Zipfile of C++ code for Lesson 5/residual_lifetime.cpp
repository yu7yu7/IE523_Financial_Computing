// Illustration of the Residual Lifetime Paradox involving 
// Exponential Distributions.
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

double get_uniform()
{
    // http://www.cplusplus.com/reference/random/uniform_real_distribution/operator()/
    std::uniform_real_distribution <double> distribution(0.0,1.0);
    double number = distribution(generator);
    return (number);
}

double get_exp(double lambda)
{
    // http://www.cplusplus.com/reference/random/exponential_distribution/
    std::exponential_distribution<double> distribution(lambda);
    double number = distribution(generator);
    return (number);
}

int main (int argc, char* argv[])
{
	int no_of_trials;
	double max_interval_length, time_traveler_arrives, lambda;
	double prev_event_time, next_event_time;
	double sum_of_sampled_intervals;
	
	sscanf (argv[1], "%d", &no_of_trials);
	sscanf (argv[2], "%lf", &max_interval_length);
	sscanf (argv[3], "%lf", &lambda);
	
	cout << "Illustration of the Residual Lifetime Paradox..." << endl;
	cout << "Number of Trials        = " << no_of_trials << endl;
	cout << "Maximum Interval Length = " << max_interval_length << endl;
	
	sum_of_sampled_intervals = 0.0;
	for (int i = 1; i <= no_of_trials; i++) {
		time_traveler_arrives = get_uniform() * max_interval_length;
		/* The traveller arrives at a time instant uniformly distributed
		 over the max_interval_length */
		
		prev_event_time = next_event_time = 0.0;
		/* run the clock till we get to the instant when the traveller arrived */
		while ((prev_event_time < time_traveler_arrives) &&
			   (next_event_time <= time_traveler_arrives)) {
			prev_event_time = next_event_time;
			next_event_time = next_event_time + get_exp(lambda);
		}
		/* use the sampled interval for future data analysis */
		sum_of_sampled_intervals = sum_of_sampled_intervals +
		(next_event_time - prev_event_time);
	}
	//cout << sum_of_sampled_intervals << endl;
	cout << "Average Computed by Random Sampling = " <<
	sum_of_sampled_intervals/((double) no_of_trials) << endl;
	cout << "Common Average                      = " << 1.0/lambda << endl;
}
