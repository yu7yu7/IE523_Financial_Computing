// Computing the Greeks (actually, just Delta; you can do the rest yourself)
// from single-run simulations using the theory of Infinitesimal Perturbation Analysis
// Written by Prof. R.S. Sreenivas for IE523: Financial Computing
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include "normdist.h" 
#define UNIT_STEP(x) ((x)>0?(1):(0))
using namespace std;

float risk_free_rate, strike_price, initial_stock_price, expiration_time, volatility;
int no_of_trials;


double option_price_put_black_scholes(const double& S,      // spot price
									  const double& K,      // Strike (exercise) price,
									  const double& r,      // interest rate
									  const double& sigma,  // volatility
									  const double& time){
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+r*time)/(sigma*time_sqrt) + 0.5*sigma*time_sqrt;
    double d2 = d1-(sigma*time_sqrt);
    return K*exp(-r*time)*N(-d2) - S*N(-d1);
};

double option_price_call_black_scholes(const double& S,       // spot (underlying) price
									   const double& K,       // strike (exercise) price,
									   const double& r,       // interest rate
									   const double& sigma,   // volatility 
									   const double& time) {  // time to maturity 
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

double option_price_delta_call_black_scholes(const double& S,     // spot price
											 const double& K,     // Strike (exercise) price,
											 const double& r,     // interest rate
											 const double& sigma, // volatility
											 const double& time){  // time to maturity
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+r*time)/(sigma*time_sqrt) + 0.5*sigma*time_sqrt; 
    double delta = N(d1);
    return delta;
};

double option_price_delta_put_black_scholes(const double& S, // spot price
											const double& K, // Strike (exercise) price,
											const double& r,  // interest rate
											const double& sigma,
											const double& time) {
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+r*time)/(sigma*time_sqrt) + 0.5*sigma*time_sqrt; 
    double delta = -N(-d1);
    return delta;
}

float max(float a, float b) {
	return (b < a )? a:b;
}

float get_uniform()
{
	return (((float) random())/(pow(2.0, 31.0)-1.0));
}

int main (int argc, char* argv[])
{	
	sscanf (argv[1], "%f", &expiration_time);
	sscanf (argv[2], "%f", &risk_free_rate);
	sscanf (argv[3], "%f", &volatility);
	sscanf (argv[4], "%f", &initial_stock_price);
	sscanf (argv[5], "%f", &strike_price);
	sscanf (argv[6], "%d", &no_of_trials);
	
	float R = (risk_free_rate - 0.5*pow(volatility,2))*expiration_time;
	float SD = volatility*sqrt(expiration_time); 
	
	cout << "--------------------------------" << endl;
	cout << "Greeks and European Options Pricing via Monte Carlo Simulation" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility*100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "Number of Trials = " << no_of_trials << endl;
	cout << "--------------------------------" << endl;
	
	// The greeks by automatic differentiation of the simulation code
	float put_delta = 0.0;
	float call_delta = 0.0; 
	float derivative_of_S1_wrt_S = 0.0;
	float derivative_of_S2_wrt_S = 0.0;
	float derivative_of_S3_wrt_S = 0.0;
	float derivative_of_S4_wrt_S = 0.0;
	
	float put_option_price = 0.0;
	float call_option_price = 0.0;
	for (int i = 0; i < no_of_trials; i++) {
		// generate unit-normals using Box-Muller Transform
		float x = get_uniform();
		float y = get_uniform();
		float a =  sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
		float b =  sqrt(-2.0*log(x)) * sin(6.283185307999998*y);
		
		// generate one humungous sample path for the underlying 
		// stock price by simulation... we could have cut-up the 
		// expiration_time into many smaller segements and produced
		// a sample path... but it would be statistically similar to
		// what we are doing below.  Anyways we only need the value at
		// expiration time for the price of an European Option.
		float S1 = initial_stock_price * exp(R + SD*a);
		derivative_of_S1_wrt_S = exp(R + SD*a);
		float S2 = initial_stock_price * exp(R - SD*a);
		derivative_of_S2_wrt_S = exp(R - SD*a);
		float S3 = initial_stock_price * exp(R + SD*b);
		derivative_of_S3_wrt_S = exp(R + SD*b);
		float S4 = initial_stock_price * exp(R - SD*b);
		derivative_of_S4_wrt_S = exp(R - SD*b);
	
		call_option_price += (max(0.0, S1 - strike_price) + 
							  max(0.0, S2 - strike_price) + 
							  max(0.0, S3 - strike_price) + 
							  max(0.0, S4 - strike_price))/4.0;
		call_delta += ((UNIT_STEP(S1-strike_price)*derivative_of_S1_wrt_S) + 
					   (UNIT_STEP(S2-strike_price)*derivative_of_S2_wrt_S) + 
					   (UNIT_STEP(S3-strike_price)*derivative_of_S3_wrt_S) + 
					   (UNIT_STEP(S4-strike_price)*derivative_of_S4_wrt_S))/4.0;
		put_option_price += (max(0.0, strike_price - S1) + 
							 max(0.0, strike_price - S2) + 
							 max(0.0, strike_price - S3) + 
							 max(0.0, strike_price - S4))/4.0;
		put_delta -= ((UNIT_STEP(strike_price - S1)*derivative_of_S1_wrt_S) + 
					  (UNIT_STEP(strike_price - S2)*derivative_of_S2_wrt_S) + 
					  (UNIT_STEP(strike_price - S3)*derivative_of_S3_wrt_S) +
					  (UNIT_STEP(strike_price - S4)*derivative_of_S4_wrt_S))/4.0;
	}
	call_option_price = exp(-risk_free_rate*expiration_time)*(call_option_price/((float) no_of_trials));
	put_option_price = exp(-risk_free_rate*expiration_time)*(put_option_price/((float) no_of_trials));
	call_delta = exp(-risk_free_rate*expiration_time)*(call_delta/((float) no_of_trials));
	put_delta = exp(-risk_free_rate*expiration_time)*(put_delta/((float) no_of_trials));
		
	cout << "The average Call Price is " << call_option_price << endl;
	cout << "Put Price according to Black-Scholes = " << 
	option_price_call_black_scholes(initial_stock_price, strike_price, risk_free_rate, 
								   volatility, expiration_time) << endl;
	cout << "Average Call Delta = " << call_delta << endl;
	cout << "Call Delta according to Black-Scholes = " <<
	option_price_delta_call_black_scholes(initial_stock_price, strike_price, risk_free_rate, 
										  volatility, expiration_time) << endl;
	cout << "The average Put Price is " << put_option_price << endl;
	cout << "Call Price according to Black-Scholes = " << 
	option_price_put_black_scholes(initial_stock_price, strike_price, risk_free_rate, 
									volatility, expiration_time) << endl;
	cout << "Average Put Delta = " << put_delta << endl;
	cout << "Put Delta according to Black-Scholes = " <<
	option_price_delta_put_black_scholes(initial_stock_price, strike_price, risk_free_rate, 
										  volatility, expiration_time) << endl;
	cout << "--------------------------------" << endl;
}	
			
		
	
