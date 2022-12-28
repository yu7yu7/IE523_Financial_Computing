#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include "normdist.h"
#include <map>
#include <random>
using namespace std;


double risk_free_rate, strike_price;
double initial_stock_price, expiration_time, volatility, barrier_price;
int no_of_trials, no_of_divisions;


unsigned seed = (unsigned) std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator(seed);

double get_uniform()
{
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    double number = distribution(generator);
    return (number);
}

float max(float a, float b) {
    return (b < a )? a:b;
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
                                       const double& time)      // time to maturity
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

double Pc(double ST)
{
    if ((ST <= barrier_price)||(initial_stock_price <= barrier_price)) return 1.0;
    else return exp( -(2*log(initial_stock_price/barrier_price)*log(ST/barrier_price) / (pow(volatility,2) * (double)expiration_time)));
}


int main (int argc, char* argv[])
{
    sscanf (argv[1], "%lf", &expiration_time);
    sscanf (argv[2], "%lf", &risk_free_rate);
    sscanf (argv[3], "%lf", &volatility);
    sscanf (argv[4], "%lf", &initial_stock_price);
    sscanf (argv[5], "%lf", &strike_price);
    sscanf (argv[6], "%d", &no_of_trials);
    sscanf (argv[7], "%d", &no_of_divisions);
    sscanf (argv[8], "%lf", &barrier_price);
    
    cout << "--------------------------------------" << endl;
    cout << "European Down-and-Out Continuous Barrier Options Pricing via Monte-Carlo Simulation" << endl;
    cout << "Expiration Time (Years) = " << expiration_time << endl;
    cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
    cout << "Volatility (%age of stock value) = " << volatility*100 << endl;
    cout << "Initial Stock Price = " << initial_stock_price << endl;
    cout << "Strike Price = " << strike_price << endl;
    cout << "Barrier Price = " << barrier_price << endl;
    cout << "Number of Trials = " << no_of_trials << endl;
    cout << "Number of Divisions = " << no_of_divisions << endl;
    cout << "--------------------------------------" << endl;
    
    double delta_T = expiration_time/((double) no_of_divisions);
    double R = (risk_free_rate - 0.5*pow(volatility,2))*delta_T;
    double SD = volatility*sqrt(delta_T);
    
    double put_option_price = 0.0;
    double call_option_price = 0.0;
    double adj_call_option_price = 0.0;
    double adj_put_option_price = 0.0;
    
    for (int i = 0; i < no_of_trials; i++) {

        double S1 = initial_stock_price;
        double S2 = initial_stock_price;
        double S3 = initial_stock_price;
        double S4 = initial_stock_price;
        
        for (int j = 0; j < no_of_divisions; j++){
            double x = get_uniform();
            double y = get_uniform();
            double a =  sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
            double b =  sqrt(-2.0*log(x)) * sin(6.283185307999998*y);
            
            if (S1 == 0.0) S1 = 0.0;
            else if (S1 * exp(R + SD*a) <= barrier_price) S1 = 0.0;
            else S1 = S1 * exp(R + SD*a);
            
            if (S2 == 0.0) S2 = 0.0;
            else if (S2 * exp(R - SD*a) <= barrier_price) S2 = 0.0;
            else S2 = S2 * exp(R - SD*a);
            
            if (S3 == 0.0) S3 = 0.0;
            else if (S3 * exp(R + SD*b) <= barrier_price) S3 = 0.0;
            else S3 = S3 * exp(R + SD*b);
            
            if (S4 == 0.0) S4 = 0.0;
            else if (S4 * exp(R - SD*b) <= barrier_price) S4 = 0.0;
            else S4 = S4 * exp(R - SD*b);
            
        }
        
        call_option_price += (max(0.0, S1 - strike_price) +
                              max(0.0, S2 - strike_price) +
                              max(0.0, S3 - strike_price) +
                              max(0.0, S4 - strike_price))/4.0;

        put_option_price += (max(0.0, strike_price - S1) +
                             max(0.0, strike_price - S2) +
                             max(0.0, strike_price - S3) +
                             max(0.0, strike_price - S4))/4.0;
        
        
        adj_call_option_price += ((max(0.0, S1 - strike_price) * (1.0-Pc(S1))) +
                              (max(0.0, S2 - strike_price) * (1.0-Pc(S2))) +
                              (max(0.0, S3 - strike_price) * (1.0-Pc(S3))) +
                              (max(0.0, S4 - strike_price) * (1.0-Pc(S4)))) / 4.0;

        adj_put_option_price += ((max(0.0, strike_price - S1) * (1 - Pc(S1))) +
                             (max(0.0, strike_price - S2) * (1 - Pc(S2))) +
                             (max(0.0, strike_price - S3) * (1 - Pc(S3))) +
                             (max(0.0, strike_price - S4) * (1 - Pc(S4)))) / 4.0;


    }
    call_option_price = exp(-risk_free_rate*expiration_time)*(call_option_price/((double) no_of_trials));
    put_option_price = exp(-risk_free_rate*expiration_time)*(put_option_price/((double) no_of_trials));
    adj_call_option_price = exp(-risk_free_rate*expiration_time)*(adj_call_option_price/((double) no_of_trials));
    adj_put_option_price = exp(-risk_free_rate*expiration_time)*(adj_put_option_price/((double) no_of_trials));
    
    cout << "--------------------------------------" << endl;
    cout << "The average Call Price by explicit simulation  = " << call_option_price << endl;
    cout << "The Call price using the (1-p)-adjustment term = " << adj_call_option_price << endl;
    cout << "Theoretical Call Price                         = " <<closed_form_down_and_out_european_call_option() << endl;
    cout << "--------------------------------------" << endl;
    cout << "The average Put Price by explicit simulation  = " << put_option_price << endl;
    cout << "The Put price using the (1-p)-adjustment term = " << adj_put_option_price << endl;
    cout << "Theoretical Put Price                         = " <<closed_form_down_and_out_european_put_option() << endl;
    cout << "--------------------------------------" << endl;

}
