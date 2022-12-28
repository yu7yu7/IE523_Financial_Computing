//
//  main.cpp
//  Gurobi European via LP
//
//  Created by Ramavarapu Sreenivas on 8/5/22.
//  Pricing an European Option using a Replicating Portfolio, which in turn
//  is found using Linear Programming (i.e. we are going to use the Gurobi
//  API).
//
//  There is a lot that you have figure out on your own in this code. There
//  is a lot of "small details" that *you* have to work-through and learn.
//  This will help you learn about a powerful optimiztion package -- Gurobi
//
//  The work of Edirisinghe, Naik and Uppal (paper can be found on Compass)
//  shows that this method can be easily extended/modified to find the price
//  of Options where there are transaction costs.
//
//  Written by Prof. Sreenivas for IE523: Financial Computation
//  It is relatively straightforward to modify this code to include transactions costs
//  as suggested in the paper.  You will have to read Lesson 6 of my notes to follow the
//  logic behind how the various constraints are ennunciated in this program.

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include "normdist.h"          // this defines the normal distribution from Odegaard's files
#include <sstream>
#include "gurobi_c++.h"

using namespace std;

string itos(int i)
{
    stringstream s;
    s << i;
    return s.str();
}


// Declaring some global variables to make my life easy
// You should look at cleaning a lot of this up...

double up_factor, risk_free_rate, strike_price, R;
double initial_stock_price, expiration_time, volatility;
int no_of_divisions;

GRBEnv call_env = GRBEnv();
GRBModel call_model = GRBModel(call_env);

GRBEnv put_env = GRBEnv();
GRBModel put_model = GRBModel(put_env);

double max(double a, double b) {
    return (b < a )? a:b;
}

// Odegaard's Black-Scholes Put C++ Code
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

// Odegaard's Black-Scholes Call C++ Code
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

// Odegaard's "Curve Fitting" of a curve into the unit normal CDF
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

// This recursive routine writes the relevant constraints for the LP
// formulation for the Call Option (Put Option is similar, mutatis-mutandis).
//
// We do not need "x" (i.e. #shares of the underlying) and "y" ($'s in "bank") for the final
// stage of the Binomial Lattice - Why?  Because, we are at expiration, and there is nothing
// afterwards.
//
// The number of "non-terminal" nodes in a Binomial Lattice with "no_of_divisions=n" is essentially
// 1 + 2 + ... + n (draw the Lattice and and check it out for yourself).  That is, there are n(n+1)/2
// "interior/non-terminal" nodes.  If (k,i) is a non-terminal node in the lattice, then it has an "x_i"
// (#shares at (k,i)) and "y_i" ($'s placed in a "bank" at (k,i)).  Therefore, we have a total of
// n*(n+1)-many (x_i, y_i) variables for the lattice.  The Gurobi API requires the row-constraint to be
// #variables, which is why all "row" arrays have a size of "(no_of_divisions*(no_of_divisions+1)"
// in the code listed below. Y

void create_LP_for_european_call_option(int k, int i, GRBVar *x)
{
    if (k == no_of_divisions-1)
    {
        // we are one-division away from expiration, so we play out both possibilities -- uptick and downtick --
        // which will end-up in us being at expiration (i.e. no "x" and "y" necessary when this happens).  That is,
        // the present "x" and "y" should be sufficiently large to make sure we are in-the-money irrespective of
        // what happens (i.e. irrespective of uptick or downtick).
        //
        // creating the constraint row & initializing with all zeros
        double row1[no_of_divisions*(no_of_divisions+1)], row2[no_of_divisions*(no_of_divisions+1)];
        for (int j = 0; j < no_of_divisions*(no_of_divisions+1); j++)
        {
            row1[j] = 0.0;
            row2[j] = 0.0;
        }
        
        // make relevant non-zero entries on this row
        // the x- and y-value's location on the constraint-row for (k,i) is
        // k^2 + k + (i+k)and k^2 + k + (i+k) + 1...
        row1[(k*k) + k + (i+k)] = initial_stock_price*pow(up_factor, ((float) i+1));
        row1[(k*k) + k + (i+k) + 1] = -R;
        row2[(k*k) + k + (i+k)] = initial_stock_price*pow(up_factor, ((float) i-1));
        row2[(k*k) + k + (i+k) + 1] = -R;
        
        GRBLinExpr expr1 = NULL;
        GRBLinExpr expr2 = NULL;
        for (int j = 0; j < no_of_divisions*(no_of_divisions+1); j++) {
            if (row1[j] != 0)
                expr1 += row1[j]*x[j];
            if (row2[j] != 0)
                expr2 += row2[j]*x[j];
        }
        call_model.addConstr(expr1 >= max(0.0, (initial_stock_price*pow(up_factor, ((float) i+1))) - strike_price));
        call_model.addConstr(expr2 >= max(0.0, (initial_stock_price*pow(up_factor, ((float) i-1))) - strike_price));
    }
    else
    {
        // creating two constraint rows & initializing with all zeros
        double row1[no_of_divisions*(no_of_divisions+1)], row2[no_of_divisions*(no_of_divisions+1)];
        for (int j = 0; j < no_of_divisions*(no_of_divisions+1); j++)
        {
            row1[j] = 0.0;
            row2[j] = 0.0;
        }
        
        // make relevant non-zero entries on this row
        row1[(k*k) + k + (i+k)] = initial_stock_price*pow(up_factor, ((float) i+1));
        row1[(k*k) + k + (i+k) + 1] = -R;
        row2[(k*k) + k + (i+k)] = initial_stock_price*pow(up_factor, ((float) i-1));
        row2[(k*k) + k + (i+k) + 1] = -R;
        
        row1[((k+1)*(k+1)) + (k+1) + (i+1+k+1)] = -initial_stock_price*pow(up_factor, ((float) i+1));
        row1[((k+1)*(k+1)) + (k+1) + (i+1+k+1) + 1] = 1;
        row2[((k+1)*(k+1)) + (k+1) + (i-1+k+1)] = -initial_stock_price*pow(up_factor, ((float) i-1));
        row2[((k+1)*(k+1)) + (k+1) + (i-1+k+1) + 1] = 1;
        
        GRBLinExpr expr1 = NULL;
        GRBLinExpr expr2 = NULL;
        for (int j = 0; j < no_of_divisions*(no_of_divisions+1); j++) {
            if (row1[j] != 0)
                expr1 += row1[j]*x[j];
            if (row2[j] != 0)
                expr2 += row2[j]*x[j];
        }
        call_model.addConstr(expr1 >= 0.0, "C"+itos(k)+itos(i));
        call_model.addConstr(expr2 >= 0.0, "D"+itos(k)+itos(i));
        
        create_LP_for_european_call_option(k+1, i+1, x);
        create_LP_for_european_call_option(k+1, i-1, x);
    }
}


// This recursive routine writes the relevant constraints for the LP
// formulation for the Put Option
void create_LP_for_european_put_option(int k, int i, GRBVar *x)
{
    if (k == no_of_divisions-1)
    {
        // creating the constraint row & initializing with all zeros
        double row1[no_of_divisions*(no_of_divisions+1)], row2[no_of_divisions*(no_of_divisions+1)];
        for (int j = 0; j < no_of_divisions*(no_of_divisions+1); j++)
        {
            row1[j] = 0.0;
            row2[j] = 0.0;
        }
        
        // make relevant non-zero entries on this row
        // the x- and y-value's location on the constraint-row for (k,i) is
        // k^2 + k + (i+k)and k^2 + k + (i+k) + 1...
        row1[(k*k) + k + (i+k)] = -initial_stock_price*pow(up_factor, ((float) i+1));
        row1[(k*k) + k + (i+k) + 1] = R;
        row2[(k*k) + k + (i+k)] = -initial_stock_price*pow(up_factor, ((float) i-1));
        row2[(k*k) + k + (i+k) + 1] = R;
        
        GRBLinExpr expr1 = NULL;
        GRBLinExpr expr2 = NULL;
        for (int j = 0; j < no_of_divisions*(no_of_divisions+1); j++) {
            if (row1[j] != 0)
                expr1 += row1[j]*x[j];
            if (row2[j] != 0)
                expr2 += row2[j]*x[j];
        }
        put_model.addConstr(expr1 >= max(0.0, (strike_price - initial_stock_price*pow(up_factor, ((float) i+1)))));
        put_model.addConstr(expr2 >= max(0.0, (strike_price - initial_stock_price*pow(up_factor, ((float) i-1)))));
    }
    else
    {
        // creating two constraint rows & initializing with all zeros
        double row1[no_of_divisions*(no_of_divisions+1)], row2[no_of_divisions*(no_of_divisions+1)];
        for (int j = 0; j < no_of_divisions*(no_of_divisions+1); j++)
        {
            row1[j] = 0.0;
            row2[j] = 0.0;
        }
        
        // make relevant non-zero entries on this row
        row1[(k*k) + k + (i+k)] = -initial_stock_price*pow(up_factor, ((float) i+1));
        row1[(k*k) + k + (i+k) + 1] = R;
        row2[(k*k) + k + (i+k)] = -initial_stock_price*pow(up_factor, ((float) i-1));
        row2[(k*k) + k + (i+k) + 1] = R;
        
        row1[((k+1)*(k+1)) + (k+1) + (i+1+k+1)] = initial_stock_price*pow(up_factor, ((float) i+1));
        row1[((k+1)*(k+1)) + (k+1) + (i+1+k+1) + 1] = -1;
        row2[((k+1)*(k+1)) + (k+1) + (i-1+k+1)] = initial_stock_price*pow(up_factor, ((float) i-1));
        row2[((k+1)*(k+1)) + (k+1) + (i-1+k+1) + 1] = -1;
        
        GRBLinExpr expr1 = NULL;
        GRBLinExpr expr2 = NULL;
        for (int j = 0; j < no_of_divisions*(no_of_divisions+1); j++) {
            if (row1[j] != 0)
                expr1 += row1[j]*x[j];
            if (row2[j] != 0)
                expr2 += row2[j]*x[j];
        }
        put_model.addConstr(expr1 >= 0.0, "E"+itos(k)+itos(i));
        put_model.addConstr(expr2 >= 0.0, "F"+itos(k)+itos(i));
        
        create_LP_for_european_put_option(k+1,i+1,x);
        create_LP_for_european_put_option(k+1,i-1,x);
    }
}

void set_up_and_solve_the_LP_for_the_call_option()
{
    // get everything started; the number of variables equals 2 * the number of
    // vertices upto the last-but-one layer before expiry
    GRBVar *x = new GRBVar[no_of_divisions*(no_of_divisions+1)];
    for (int i = 0; i < no_of_divisions*(no_of_divisions+1); i++)
        x[i] = call_model.addVar(0.0, GRB_INFINITY, 0, GRB_CONTINUOUS, "C"+itos(i));
    
    // This keeps the message reporting of Gurobi to a minimum
    call_model.set(GRB_IntParam_OutputFlag, 0);
    
    // set the constraints in the LP for a call option
    create_LP_for_european_call_option(0, 0, x);
    
    // set the objective function
    GRBLinExpr call_obj = initial_stock_price*x[0] - x[1];
    call_model.setObjective(call_obj);
    
    //call_model.write("call.lp");
    
    // solve the LP
    call_model.optimize();
    
    {
        int optimstatus = call_model.get(GRB_IntAttr_Status);
        if (optimstatus == GRB_OPTIMAL)
            cout << "Call Price according the LP formulation = " << call_model.get(GRB_DoubleAttr_ObjVal) << endl;
        else
            cout << "Gurobi ran into problems;  GRB_IntAttr_Status = " << optimstatus << endl;
    }
}

void set_up_and_solve_the_LP_for_the_put_option()
{
    // get everything started; the number of variables equals 2 * the number of
    // vertices upto the last-but-one layer before expiry
    GRBVar *x = new GRBVar[no_of_divisions*(no_of_divisions+1)];
    for (int i = 0; i < no_of_divisions*(no_of_divisions+1); i++)
        x[i] = put_model.addVar(0.0, GRB_INFINITY, 0, GRB_CONTINUOUS, "P"+itos(i));
    
    // This keeps the message reporting of Gurobi to a minimum
    put_model.set(GRB_IntParam_OutputFlag, 0);
    
    // set the constraints in the LP for a call option
    create_LP_for_european_put_option(0, 0, x);
    
    // set the objective function
    GRBLinExpr put_obj = x[1] - initial_stock_price*x[0];
    put_model.setObjective(put_obj);
    
    //put_model.write("put.lp");
    
    // solve the LP
    put_model.optimize();
    
    {
        int optimstatus = put_model.get(GRB_IntAttr_Status);
        if (optimstatus == GRB_OPTIMAL)
            cout << "Put Price according the LP formulation = " << put_model.get(GRB_DoubleAttr_ObjVal) << endl;
        else
            cout << "Gurobi ran into problems;  GRB_IntAttr_Status = " << optimstatus << endl;
    }
}

int main (int argc, char* argv[])
{
    sscanf (argv[1], "%lf", &expiration_time);
    sscanf (argv[2], "%d", &no_of_divisions);
    sscanf (argv[3], "%lf", &risk_free_rate);
    sscanf (argv[4], "%lf", &volatility);
    sscanf (argv[5], "%lf", &initial_stock_price);
    sscanf (argv[6], "%lf", &strike_price);
    
    up_factor = exp(volatility*sqrt(expiration_time/((double) no_of_divisions)));
    R = exp(risk_free_rate*expiration_time/((double) no_of_divisions));
    
    cout << "--------------------------------------" << endl;
    cout << "European Option Pricing via Linear Programming" << endl;
    cout << "Expiration Time (Years) = " << expiration_time << endl;
    cout << "Number of Divisions = " << no_of_divisions << endl;
    cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
    cout << "Volatility (%age of stock value) = " << volatility*100 << endl;
    cout << "Initial Stock Price = " << initial_stock_price << endl;
    cout << "Strike Price = " << strike_price << endl;
    cout << "R = " << R << endl;
    cout << "Up-Factor = " << up_factor << endl;
    cout << "--------------------------------------" << endl;

    set_up_and_solve_the_LP_for_the_call_option();
    cout << "Call Price according to Black-Scholes = " <<
    option_price_call_black_scholes(initial_stock_price, strike_price, risk_free_rate,
                                    volatility, expiration_time) << endl;
    cout << "--------------------------------------" << endl;
    
    set_up_and_solve_the_LP_for_the_put_option();
    cout << "Put Price according to Black-Scholes = " <<
    option_price_put_black_scholes(initial_stock_price, strike_price, risk_free_rate,
                                    volatility, expiration_time) << endl;
    cout << "--------------------------------------" << endl;

}
 
