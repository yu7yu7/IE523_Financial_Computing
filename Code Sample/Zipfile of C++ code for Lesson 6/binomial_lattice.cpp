//
//  main.cpp
//  Binomial Lattice
//
//  Created by Ramavarapu Sreenivas on 11/4/18.
//  Copyright © 2018 Ramavarapu Sreenivas. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <random>
#include <chrono>
using namespace std;
#define PI 3.141592654

double risk_free_rate, initial_stock_price, expiration_time, volatility;
int no_of_stages;

long **memoized_n_take_i;

unsigned seed = (unsigned) std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator (seed);

// u.i.i.d. generator
double get_uniform()
{
    std::uniform_real_distribution <double> distribution(0.0,1.0);
    double number = distribution(generator);
    return (number);
}

// unit-normal i.i.d. generator
double get_gaussian()
{
    return (sqrt(-2.0*log(get_uniform()))*cos(6.283185307999998*get_uniform()));
}

int take(int n, int i)
{
    if ((n < 0) || (i < 0) || (i > n))
        return (0);
    if ((0 == i) || (n == i))
        return (1);
    if (1 == i)
        return (n);
    return (take(n-1,i-1) + take(n-1,i));
}

double prob_of_final_asset_price(double x)
{
    // parameters for the log-normal terminal asset price
    // I will need this for the log-normal-pdf (equation 4; lesson 6)
    double R = (risk_free_rate - 0.5*pow(volatility,2))*expiration_time;
    double SD = volatility*sqrt(expiration_time);
    double result = (1.0/sqrt(2*PI))*(1.0/(volatility*sqrt(expiration_time)))*(1.0/x);
    result *= exp(-(log(x)-log(initial_stock_price)-R)*(log(x)-log(initial_stock_price)-R)/(2.0*SD*SD));
    return (result);
}

int main (int argc, char* argv[])
{
    sscanf (argv[1], "%lf", &expiration_time);
    sscanf (argv[2], "%lf", &risk_free_rate);
    sscanf (argv[3], "%lf", &volatility);
    sscanf (argv[4], "%lf", &initial_stock_price);
    sscanf (argv[5], "%d", &no_of_stages);
    ofstream outf(argv[6]);
    
    memoized_n_take_i = new long*[no_of_stages+1];
    for (int i = 0; i <= no_of_stages; i++)
        memoized_n_take_i[i] = new long[no_of_stages+1];
    for (int i = 0; i <= no_of_stages; i++)
        for (int j = 0; j <= no_of_stages; j++)
            memoized_n_take_i[i][j] = -1;
    
    // parameters for the binomial lattice (I am using lattice_R here; not to be confused with R
    // from the log-normal asset price parameters
    double up_factor = exp(volatility*sqrt(expiration_time/((double) no_of_stages)));
    double lattice_R = exp(risk_free_rate*expiration_time/((float) no_of_stages));
    double uptick_prob = (lattice_R - (1/up_factor))/(up_factor-(1/up_factor));
    double downtick_prob = 1.0 - uptick_prob;
    
    cout << "--------------------------------" << endl;
    cout << "Asset-Price Distribution at Expiration (Binomial Lattice vs. Log-Normal)" << endl;
    cout << "Expiration Time (Years) = " << expiration_time << endl;
    cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
    cout << "Volatility (%age of stock value) = " << volatility*100 << endl;
    cout << "Initial Stock Price = " << initial_stock_price << endl;
    cout << "Number of Lattice Stages = " << no_of_stages << endl;
    cout << "Uptick Probability = " << uptick_prob << endl;
    cout << "--------------------------------" << endl;
    
    double price_vector[no_of_stages+1];
    double probability_vector[no_of_stages+1];
    double normalized_log_normal_probabilities[no_of_stages+1];
    double sum = 0.0;
    // The terminal-stage of the binomial lattice will have (#stages + 1) nodes, the i-th node
    // (i = 0,1,...,#stages) has an associated-price of initial_stock_price*(upfactor)^(2*i-#stages)
    // The probability-mass (n-take-i) (downtick_prob)^(n-i) (uptick_prob)^i
    for (int i = 0; i <= no_of_stages; i++)
    {
        price_vector[i] = initial_stock_price*pow(up_factor, ((double) ((2*i) - no_of_stages)));
        probability_vector[i] = take(no_of_stages,i)*pow(downtick_prob, (double) no_of_stages-i)*pow(uptick_prob, (double) i);
        normalized_log_normal_probabilities[i] = prob_of_final_asset_price(price_vector[i]);
        sum += normalized_log_normal_probabilities[i];
    }
    
    for (int i = 0; i <= no_of_stages; i++)
        outf << price_vector[i] << " " << probability_vector[i] << " " << normalized_log_normal_probabilities[i]/sum << endl;

}
