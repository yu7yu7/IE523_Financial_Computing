// Down-and-out European Barrier Option Pricing via truncated binomial lattices
// Written by Prof. Sreenivas
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include "normdist.h"
using namespace std;

float up_factor, uptick_prob, risk_free_rate, strike_price, barrier_price;
float initial_stock_price, expiration_time, volatility, R;
int no_of_divisions;

float max(float a, float b) {
	return (b < a )? a:b;
}

float european_down_and_out_call_option(int k, int i) {
	if (k == no_of_divisions) {
		if (initial_stock_price*pow(up_factor, ((float) i)) > barrier_price)
			return max(0.0, (initial_stock_price*pow(up_factor, ((float) i))) - strike_price);
		else
			return (0.0);
	}
	else {
		if (initial_stock_price*pow(up_factor, ((float) i)) > barrier_price)
			return ((uptick_prob*european_down_and_out_call_option(k+1,i+1) +
					 (1-uptick_prob)*european_down_and_out_call_option(k+1,i-1))/R);
		else
			return (0.0);
	}
}

float european_down_and_out_put_option(int k, int i) {
	if (k == no_of_divisions) {
		if (initial_stock_price*pow(up_factor, ((float) i)) > barrier_price)
			return max(0.0, (strike_price - initial_stock_price*pow(up_factor, ((float) i))));
		else
			return (0.0);
	}
	else {
		if (initial_stock_price*pow(up_factor, ((float) i)) > barrier_price)
			return ((uptick_prob*european_down_and_out_put_option(k+1,i+1) +
					 (1-uptick_prob)*european_down_and_out_put_option(k+1,i-1))/R);
		else
			return (0.0);
	}
}

double option_price_put_black_scholes(const double& S,      // spot price
									  const double& K,      // Strike (exercise) price,
									  const double& r,      // interest rate
									  const double& sigma,  // volatility
									  const double& time)
{
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+r*time)/(sigma*time_sqrt) + 0.5*sigma*time_sqrt;
    double d2 = d1-(sigma*time_sqrt);
    return K*exp(-r*time)*N(-d2) - S*N(-d1);
};

double option_price_call_black_scholes(const double& S,       // spot (underlying) price
									   const double& K,       // strike (exercise) price,
									   const double& r,       // interest rate
									   const double& sigma,   // volatility 
									   const double& time)	  // time to maturity 
{  
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+r*time)/(sigma*time_sqrt)+0.5*sigma*time_sqrt; 
    double d2 = d1-(sigma*time_sqrt);
    return S*N(d1) - K*exp(-r*time)*N(d2);
};

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

float closed_form_down_and_out_european_call_option()
{
	// I took this formula from Wilmott, Howison and Dewynne, "The Mathematics of Financial Derivatives"
	float K = (2*risk_free_rate)/(volatility*volatility);
	float A = option_price_call_black_scholes(initial_stock_price, strike_price, 
											  risk_free_rate, volatility, expiration_time);
	float B = (barrier_price*barrier_price)/initial_stock_price;
	float C = pow(initial_stock_price/barrier_price, -(K-1));
	float D = option_price_call_black_scholes(B, strike_price, risk_free_rate, volatility, expiration_time);
	return (A - D*C);
}

float closed_form_down_and_in_european_put_option()	
{
	// just making it easier by renaming the global variables locally
	float S = initial_stock_price;
	float r = risk_free_rate;
	float T = expiration_time;
	float sigma = volatility;
	float H = barrier_price;
	float X = strike_price;
	
	// Took these formulae from some online reference
	float lambda = (r+((sigma*sigma)/2))/(sigma*sigma);
	float temp = 2*lambda - 2.0;
	float x1 = (log(S/H)/(sigma*sqrt(T))) + (lambda*sigma*sqrt(T));
	float y = (log(H*H/(S*X))/(sigma*sqrt(T))) + (lambda*sigma*sqrt(T));
	float y1 = (log(H/S)/(sigma*sqrt(T))) + (lambda*sigma*sqrt(T));
	return (-S*N(-x1) + X*exp(-r*T)*N(-x1 + sigma*sqrt(T)) + 
			S*pow(H/S, 2*lambda)*(N(y)-N(y1)) -
			X*exp(-r*T)*pow(H/S, temp)*(N(y-sigma*sqrt(T))-N(y1-sigma*sqrt(T))));
}

float closed_form_down_and_out_european_put_option()
{
	float vanilla_put = option_price_put_black_scholes(initial_stock_price, strike_price, 
													   risk_free_rate, volatility, expiration_time);
	float put_down_in = closed_form_down_and_in_european_put_option();
	return (vanilla_put - put_down_in);
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
	
	up_factor = exp(volatility*sqrt(expiration_time/((float) no_of_divisions)));
	R = exp(risk_free_rate*expiration_time/((float) no_of_divisions));
	uptick_prob = (R - (1/up_factor))/(up_factor-(1/up_factor));
	
	cout << "European Down and Out Option Pricing" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Number of Divisions = " << no_of_divisions << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility*100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "Barrier Price = " << barrier_price << endl; 
	cout << "--------------------------------------" << endl;
	cout << "Up Factor = " << up_factor << endl;
	cout << "Uptick Probability = " << uptick_prob << endl;
	cout << "--------------------------------------" << endl;
	cout << "Binomial Price of an European Down and Out Call Option = " << european_down_and_out_call_option(0, 0) << endl;
	cout << "Binomial Price of an European Down and Out Put Option = " << european_down_and_out_put_option(0, 0) << endl;
	cout << "--------------------------------------" << endl;
	cout << "Price of an European Down and Out Call Option from Theory = " <<
	closed_form_down_and_out_european_call_option() << endl;
	cout << "Price of an European Down and Out Put Option from Theory = " <<
	closed_form_down_and_out_european_put_option() << endl;
	cout << "--------------------------------------" << endl;
}