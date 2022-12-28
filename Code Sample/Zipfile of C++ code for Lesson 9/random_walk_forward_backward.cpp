// Simulating a random walk and then trying to go backwards after that
// Written by Prof. Sreenivas for IE523: Financial Computing
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
using namespace std;

float risk_free_rate, initial_stock_price, expiration_time, volatility;
int no_of_trials, no_of_divisions;

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
	sscanf (argv[5], "%d", &no_of_divisions);
	ofstream outf(argv[6]); // 6th command line argument is the CSV output file name
	
	// chop expiration time into no_of_divisions many segments 
	// figure out the motion within each segment
	float delta_T = expiration_time/((float) no_of_divisions);
	float delta_R = (risk_free_rate - 0.5*pow(volatility,2))*delta_T;
	float delta_SD = volatility*sqrt(delta_T); 

	// by sharing random variables we create 4 paths 
	float current_stock_price1 = initial_stock_price;
	float current_stock_price2 = initial_stock_price;
	float current_stock_price3 = initial_stock_price;
	float current_stock_price4 = initial_stock_price;
	for (int j = 0; j < no_of_divisions; j++) {
		
		outf << expiration_time/((float) no_of_divisions) * ((float) j) << ", " <<
		current_stock_price1 << ", " << current_stock_price2 << ", " <<
		current_stock_price3 << ", " << current_stock_price4 << endl;
		
		// create the unit normal variates using the Box-Muller Transform
		float x = get_uniform();
		float y = get_uniform();
		float a =  sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
		float b =  sqrt(-2.0*log(x)) * sin(6.283185307999998*y);
		
		current_stock_price1 = current_stock_price1*exp(delta_R + delta_SD*a);
		current_stock_price2 = current_stock_price2*exp(delta_R - delta_SD*a);
		current_stock_price3 = current_stock_price3*exp(delta_R + delta_SD*b);
		current_stock_price4 = current_stock_price4*exp(delta_R - delta_SD*b);
	}
	
	// let us see what happens when we try to run the simulation backwards from 
	// the final value... it the plot time will start from 1 to T (forward-path)
	// and then from T to 2T (backward-path)
	
	for (int j = 0; j < no_of_divisions; j++) {
		
		outf << expiration_time + (expiration_time/((float) no_of_divisions) * ((float) j)) << ", " <<
		current_stock_price1 << ", " << current_stock_price2 << ", " <<
		current_stock_price3 << ", " << current_stock_price4 << endl;
		
		// create the unit normal variates using the Box-Muller Transform
		float x = get_uniform();
		float y = get_uniform();
		float a =  sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
		float b =  sqrt(-2.0*log(x)) * sin(6.283185307999998*y);
		
		current_stock_price1 = current_stock_price1*exp(-delta_R + delta_SD*a);
		current_stock_price2 = current_stock_price2*exp(-delta_R - delta_SD*a);
		current_stock_price3 = current_stock_price3*exp(-delta_R + delta_SD*b);
		current_stock_price4 = current_stock_price4*exp(-delta_R - delta_SD*b);
	}
	
	
	cout << "--------------------------------" << endl;
	cout << "Forward-Backward Sample Path Synthesis" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility*100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Number of Steps/Divisions in the sample path = " << no_of_divisions << endl;
	cout << "The CSV output filename: " << argv[6] << endl;
	cout << "--------------------------------" << endl;
}
