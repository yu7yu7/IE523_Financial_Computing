// Down-and-out European Discrete Barrier Option Pricing Code written by Prof. Sreenivas
// The terminal payoff vector that corresponds to indices of the underlying 
// stock prices that are in the black vis-a-vis the barrier have been discounted 
// appropriately using the Brownian Bridge adjustment
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include "/Users/yu-chingliao/Desktop/Others/newmat10/newmat.h"
#include "normdist.h"
using namespace std;

float up_factor, uptick_prob, risk_free_rate, strike_price;
float initial_stock_price, expiration_time, volatility, R, barrier_price;
int no_of_divisions, no_of_sampling_instants;

float max(float a, float b) {
	return (b < a )? a:b;
}


double N(const double& z) { 
	if (z > 6.0) { return 1.0; }; // this guards against overflow 
	if (z < -6.0) { return 0.0; }; 
	double b1 = 0.31938153; 
	double b2 = -0.356563782; 
	double b3 = 1.781477937; 
	double b4 = -1.821255978; 
	double b5 = 1.330274429; 
	double p = 0.2316419; 
	double c2 = 0.3989423; 
	double a=fabs(z); 
	double t = 1.0/(1.0+a*p); 
	double b = c2*exp((-z)*(z/2.0)); 
	double n = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t; 
	n = 1.0-b*n; 
	if ( z < 0.0 ) n = 1.0 - n; 
	return n; 
}; 

Matrix repeated_squaring(Matrix A, int exponent, int no_rows)
{
	if (exponent == 0) {
		IdentityMatrix I(no_rows);
		return (I);
	}
	else if (exponent%2 == 0)
		return (repeated_squaring(A*A, exponent/2, no_rows));
	else 
		return (A*repeated_squaring(A*A, (exponent-1)/2, no_rows));
}

float european_discrete_down_and_out_call_option()
{
	int threshold;
	
	// the smallest index i for which initial_stock_price * up_factor^i exceeds barrier_price
	if (log(barrier_price/initial_stock_price)/log(up_factor) < 0.0)
		threshold = log(barrier_price/initial_stock_price)/log(up_factor);
	else
		threshold = log(barrier_price/initial_stock_price)/log(up_factor) + 1.0;
	
	// create the probability matrix
	Matrix transition_probability(2*no_of_divisions-1, 2*no_of_divisions-1);
	Matrix final_matrix(2*no_of_divisions-1, 2*no_of_divisions-1);
	
	for (int i = 1; i <= 2*no_of_divisions-1; i++)
		for (int j = 1; j <= 2*no_of_divisions-1; j++)
			transition_probability(i,j) = 0.0;
	
	// boundary values of the probabilities need to be entered
	transition_probability(1,2) = 1.0;
	transition_probability(2*no_of_divisions-1,2*no_of_divisions-1) = 1.0;
	
	for (int i = 2; i <= 2*no_of_divisions-2; i++)
		for (int j = 1; j <= 2*no_of_divisions-1; j++) {
			if (j == (i-1))
				transition_probability(i,j) = 1 - uptick_prob;
			if (j == (i+1))
				transition_probability(i,j) = uptick_prob;
		}
	
	// compute pi^no_of_divisions by method of repeated squaring
	// keep in mind that we are computing Pi^(T-1) not Pi^T 
	// because the PI matrix in itself represents the 1-step state
	// transition probability... for k-step you need to multiply PI
	// with itself (k-1) times 
	final_matrix = repeated_squaring(transition_probability, no_of_divisions-1, 2*no_of_divisions-1);
	
	// value at expiration with Bardon-Adesi, Fusari and Theal adjustments applied to all valuations
	// in the black... the ones in the red are zeroed out
	Matrix V_T(2*no_of_divisions-1,1);
	Matrix mean_at_sampling_instant(no_of_sampling_instants,1);
	Matrix variance_at_sampling_instant(no_of_sampling_instants,1);
	for (int i = 1; i <= 2*no_of_divisions-1; i++) {
		if (initial_stock_price*pow(up_factor,i-no_of_divisions) > barrier_price) {
			// Brownian Bridge w(0) = initial_stock_price, w(T) = initial_stock_price*pow(up_factor, i-no_of_divisions)
			// If we have a discrete knock out barrier on date t, then we need to compute 
			// (no_of_divisions/no_of_discrete_knock_out_barriers)-many means and variances and then use the cummulative 
			// normal to figure out the probability that the stock price at the discrete sampling time was below the 
			// barrier level.
			for (int j = 1; j <= no_of_sampling_instants; j++) {
				mean_at_sampling_instant(j,1) = initial_stock_price + 
				( ((float) j)/((float) no_of_sampling_instants)*initial_stock_price*(pow(up_factor, i-no_of_divisions)-1.0) );
				variance_at_sampling_instant(j,1) = ( ((float) j)/((float) no_of_sampling_instants) )*expiration_time*
				(1.0 - ((float) j)/((float) no_of_sampling_instants));
			}
			
			float probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant = 1.0;
			for (int j = 1; j <= no_of_sampling_instants; j++) 
				probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant *=
				(1.0 - N((barrier_price - mean_at_sampling_instant(j,1))/sqrt(variance_at_sampling_instant(j,1))));
				
			V_T(i,1) = (1.0/pow(R, no_of_divisions)) * 
			max(0.0, (initial_stock_price*pow(up_factor,i-no_of_divisions)) - strike_price) * 
			probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant;
		}
		else
			V_T(i,1) = 0.0;
	}
	
	float final_result = 0;
	for (int i = 1; i < 2*no_of_divisions-1; i++)
		final_result += final_matrix(no_of_divisions,i)*V_T(i,1);
	
	return (final_result);
}

float european_discrete_down_and_out_put_option()
{
	int threshold;
	
	// the smallest index i for which initial_stock_price * up_factor^i exceeds barrier_price
	if (log(barrier_price/initial_stock_price)/log(up_factor) < 0.0)
		threshold = log(barrier_price/initial_stock_price)/log(up_factor);
	else
		threshold = log(barrier_price/initial_stock_price)/log(up_factor) + 1.0;
	
	// create the probability matrix
	Matrix transition_probability(2*no_of_divisions-1, 2*no_of_divisions-1);
	Matrix final_matrix(2*no_of_divisions-1, 2*no_of_divisions-1);
	
	for (int i = 1; i <= 2*no_of_divisions-1; i++)
		for (int j = 1; j <= 2*no_of_divisions-1; j++)
			transition_probability(i,j) = 0.0;
	
	// boundary values of the probabilities need to be entered
	transition_probability(1,2) = 1.0;
	transition_probability(2*no_of_divisions-1,2*no_of_divisions-1) = 1.0;
	
	for (int i = 2; i <= 2*no_of_divisions-2; i++)
		for (int j = 1; j <= 2*no_of_divisions-1; j++) {
			if (j == (i-1))
				transition_probability(i,j) = 1 - uptick_prob;
			if (j == (i+1))
				transition_probability(i,j) = uptick_prob;
		}
	
	// compute pi^no_of_divisions by method of repeated squaring
	// keep in mind that we are computing Pi^(T-1) not Pi^T 
	// because the PI matrix in itself represents the 1-step state
	// transition probability... for k-step you need to multiply PI
	// with itself (k-1) times 
	final_matrix = repeated_squaring(transition_probability, no_of_divisions-1, 2*no_of_divisions-1);
	
	// value at expiration with Bardon-Adesi, Fusari and Theal adjustments applied to all valuations
	// in the black... the ones in the red are zeroed out
	Matrix V_T(2*no_of_divisions-1,1);
	Matrix mean_at_sampling_instant(no_of_sampling_instants,1);
	Matrix variance_at_sampling_instant(no_of_sampling_instants,1);
	for (int i = 1; i <= 2*no_of_divisions-1; i++) {
		if (initial_stock_price*pow(up_factor,i-no_of_divisions) > barrier_price) {
			// Brownian Bridge w(0) = initial_stock_price, w(T) = initial_stock_price*pow(up_factor, i-no_of_divisions)
			// If we have a discrete knock out barrier on date t, then we need to compute 
			// (no_of_divisions/no_of_discrete_knock_out_barriers)-many means and variances and then use the cummulative 
			// normal to figure out the probability that the stock price at the discrete sampling time was below the 
			// barrier level.
			for (int j = 1; j <= no_of_sampling_instants; j++) {
				mean_at_sampling_instant(j,1) = initial_stock_price + 
				( ((float) j)/((float) no_of_sampling_instants)*initial_stock_price*(pow(up_factor, i-no_of_divisions)-1.0) );
				variance_at_sampling_instant(j,1) = ( ((float) j)/((float) no_of_sampling_instants) )*expiration_time*
				(1.0 - ((float) j)/((float) no_of_sampling_instants));
			}
			
			float probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant = 1.0;
			for (int j = 1; j <= no_of_sampling_instants; j++) 
				probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant *=
				(1.0 - N((barrier_price - mean_at_sampling_instant(j,1))/sqrt(variance_at_sampling_instant(j,1))));
			
			V_T(i,1) = (1.0/pow(R, no_of_divisions)) * 
			max(0.0, (strike_price - initial_stock_price*pow(up_factor,i-no_of_divisions))) * 
			probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant;
		}
		else
			V_T(i,1) = 0.0;
	}
	
	float final_result = 0;
	for (int i = 1; i < 2*no_of_divisions-1; i++)
		final_result += final_matrix(no_of_divisions,i)*V_T(i,1);
	
	return (final_result);
}


int main (int argc, char* argv[])
{
	
	sscanf (argv[1], "%f", &expiration_time);
	sscanf (argv[2], "%d", &no_of_divisions);
	sscanf (argv[3], "%f", &risk_free_rate);
	sscanf (argv[4], "%f", &volatility);
	sscanf (argv[5], "%f", &initial_stock_price);
	sscanf (argv[6], "%f", &strike_price);
	sscanf (argv[7], "%f", &barrier_price);
	sscanf (argv[8], "%d", &no_of_sampling_instants);
	
	up_factor = exp(volatility*sqrt(expiration_time/((float) no_of_divisions)));
	R = exp(risk_free_rate*expiration_time/((float) no_of_divisions));
	uptick_prob = (R - (1/up_factor))/(up_factor-(1/up_factor));
	
	cout << "European Down and Out Discrete Barrier Option Pricing" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Number of Divisions = " << no_of_divisions << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility*100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "Barrier Price = " << barrier_price << endl; 
	cout << "Number of (uniformly spaced) Discrete Barrier-Samples = " << no_of_sampling_instants << endl;
	cout << "--------------------------------------" << endl;
	cout << "Up Factor = " << up_factor << endl;
	cout << "Uptick Probability = " << uptick_prob << endl;
	cout << "--------------------------------------" << endl;
	cout << "Price of an European Down and Out Discrete Barrier Call Option = " << european_discrete_down_and_out_call_option() << endl;
	cout << "Price of an European Down and Out Discrete Barrier Put Option = " << european_discrete_down_and_out_put_option() << endl;
	cout << "--------------------------------------" << endl;
}
