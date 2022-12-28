// Pricing an american put option using Longstaff and Schwartz's Least Squares Monte Carlo 
// method.    
// Written by Prof. Sreenivas for IE523: Financial Computing
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include "newmat.h"

float risk_free_rate, strike_price, initial_stock_price, expiration_time, volatility;
int no_of_trials, no_of_divisions, no_of_simulations;

float max(float a, float b) {
	return (b < a )? a:b;
}

int min(int a, int b) {
	return (b < a )? b:a;
}

float get_uniform()
{
	return (((float) random())/(pow(2.0, 31.0)-1.0));
}

float get_gaussian()
{
	return (sqrt(-2.0*log(get_uniform()))*cos(2*3.141592654*get_uniform()));
}

Matrix polynomial_regression(Matrix Independent_Variables, Matrix Dependent_Variable, int order, int no_of_observations)
{
	Matrix X(no_of_observations, order);
	Matrix Y(no_of_observations, 1);
	
	for (int i = 1; i <= no_of_observations; i++)
		Y(i,1) = Dependent_Variable(i,1);
	
	for (int j = 1; j <= order; j++) 
		for (int i = 1; i <= no_of_observations; i++) 
			X(i,j) = pow(Independent_Variables(i,1), j-1);
	
	// return inv(XT*X)*XT*Y
	Matrix X_transpose_times_X(order, order);
	X_transpose_times_X = X.t()*X;
	return (X_transpose_times_X.i() * X.t() * Y);
}


int main (int argc, char* argv[])
{
	
	sscanf (argv[1], "%f", &expiration_time);
	sscanf (argv[2], "%f", &risk_free_rate);
	sscanf (argv[3], "%f", &volatility);
	sscanf (argv[4], "%f", &initial_stock_price);
	sscanf (argv[5], "%f", &strike_price);
	sscanf (argv[6], "%d", &no_of_divisions);
	sscanf (argv[7], "%d", &no_of_simulations);
	
	if (no_of_divisions > 1000) {
		cout << "The program is currently coded for no_of_divisions <= 1000" << endl;
		cout << "Exiting..." << endl;
		return (0);
	}
	
	float delta_T = expiration_time/((float) no_of_divisions);
	float delta_R = (risk_free_rate - 0.5*pow(volatility,2))*delta_T;
	float delta_SD = volatility*sqrt(delta_T); 
	float R = exp(risk_free_rate*expiration_time/((float) no_of_divisions));
	
	cout << "--------------------------------" << endl;
	cout << "American Put Option Pricing using Longstaff and Schwartz's Least Squares Monte Carlo Simulation" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility*100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "Number of Simulations = " << no_of_simulations << endl;
	cout << "Number of Divisions = " << no_of_divisions << endl;
	cout << "R = " << R << endl;
	cout << "--------------------------------" << endl;
	
	// given array size limitations, I will run batches of 200 runs if no_of_simulations 
	// exceeds 200.  
	
	float put_price = 0.0;
	for (int k = 0; k < (no_of_simulations/200 + 1); k++) {
		
		if (k != no_of_simulations/200) 
			no_of_trials = 200;
		else 
			no_of_trials = no_of_simulations%200;
		
		if (no_of_trials != 0) {
			float asset_price[200][1000];
			
			for (int i = 0; i < no_of_trials; i++)
				asset_price[i][0] = initial_stock_price;
			
			for (int i = 0; i < no_of_trials; i++) 
				for (int j = 1; j < no_of_divisions; j++) 
					asset_price[i][j] = asset_price[i][j-1]*exp(delta_R + delta_SD*get_gaussian());

			float value[no_of_trials];
			// initialize the value based on the price at final stage
			for (int i = 0; i < no_of_trials; i++) 
				value[i] = max(0.0, strike_price - asset_price[i][no_of_divisions-1]);
			
			for (int i = (no_of_divisions-1); i > 0; i--) {
				Matrix independent_variables(no_of_trials,1);
				Matrix dependent_variables(no_of_trials,1);
				int no_of_variables = 0;
				for (int j = 0; j < no_of_trials; j++) {
					if (max(0.0, strike_price - asset_price[j][i]) > 0) {
						no_of_variables++;
						independent_variables(no_of_variables, 1) = asset_price[j][i];
						dependent_variables(no_of_variables, 1) = value[j]/R;
					}
				}
				
				if (no_of_variables > 0) {
					
					// regressing the dependent_variables on the independent variables using a 4th order polynomial
					Matrix a(min(5,no_of_variables),1);
					a = polynomial_regression(independent_variables, dependent_variables, min(5,no_of_variables), no_of_variables);
					if (no_of_variables >= 5) {
						for (int j = 0; j < no_of_trials; j++) {
							if ( ((strike_price - asset_price[j][i]) > (a(1,1) + (a(2,1)*asset_price[j][i]) + (a(3,1)*pow(asset_price[j][i],2)) + (a(4,1)*pow(asset_price[j][i],3)) + 
																		(a(5,1)*pow(asset_price[j][i],4)))) && 
								( (strike_price -asset_price[j][i]) > 0.0 ) )
								value[j] = strike_price - asset_price[j][i];
							else
								value[j] = value[j]/R;
						}
					}
					else if (no_of_variables == 4) {
						for (int j = 0; j < no_of_trials; j++) {
							if ( ((strike_price - asset_price[j][i]) > (a(1,1) + (a(2,1)*asset_price[j][i]) + (a(3,1)*pow(asset_price[j][i],2)) + (a(4,1)*pow(asset_price[j][i],3)))) && 
								( (strike_price -asset_price[j][i]) > 0.0 ) )
								value[j] = strike_price - asset_price[j][i];
							else
								value[j] = value[j]/R;
						}
						
					}
					else if (no_of_variables == 3) {
						for (int j = 0; j < no_of_trials; j++) {
							if ( ((strike_price - asset_price[j][i]) > (a(1,1) + (a(2,1)*asset_price[j][i]) + (a(3,1)*pow(asset_price[j][i],2)))) &&
								( (strike_price -asset_price[j][i]) > 0.0 ) )
								value[j] = strike_price - asset_price[j][i];
							else
								value[j] = value[j]/R;
						}	
					}
					else if (no_of_variables == 2) {
						for (int j = 0; j < no_of_trials; j++) {
							if ( ((strike_price - asset_price[j][i]) > (a(1,1) + (a(2,1)*asset_price[j][i]))) &&
								( (strike_price -asset_price[j][i]) > 0.0 ) )
								value[j] = strike_price - asset_price[j][i];
							else
								value[j] = value[j]/R;
						}
					}
					else  {
						for (int j = 0; j < no_of_trials; j++) {
							if ( ((strike_price - asset_price[j][i]) > a(1,1))  &&
								( (strike_price -asset_price[j][i]) > 0.0 ) )
								value[j] = strike_price - asset_price[j][i];
							else
								value[j] = value[j]/R;
						}	
					}
					
				}
			}
			
			float local_put_price = 0.0;
			for (int j = 0; j < no_of_trials; j++) 
				local_put_price += value[j];
			local_put_price = (local_put_price/((float) no_of_trials))/R;
			put_price += local_put_price;
		}
	}
	
	if (no_of_simulations%200 == 0) 
		cout << "American Put Price = " << put_price/((float) no_of_simulations/200)  << endl;
	else 
		cout << "Put Price = " << put_price/((float) no_of_simulations/200 + 1) << endl;
}
