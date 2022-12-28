// Multivariate Gaussian generator
// Written by Prof. Sreenivas for IE523: Financial Computing
// You will need the NEWMAT library for this to work
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <random>
#include <chrono>
#include "include.h"
#include "newmat.h"

using namespace std;

// cf http://www.cplusplus.com/reference/random/uniform_real_distribution/operator()/
// If you want to set a seed -- do it only after debug phase is completed
// otherwise errors will not be repeatable.
unsigned seed = (unsigned) std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator (seed);

ReturnMatrix Cholesky(const SymmetricMatrix& S);

double get_gaussian(double mean, double standard_deviation)
{
    std::normal_distribution<double> distribution(mean, standard_deviation);
    double number = distribution(generator);
    return (number);
}

ColumnVector multivariate_gaussian(int n, SymmetricMatrix C, ColumnVector mean)
{
    //compute the Cholesky decomposition of C
    LowerTriangularMatrix L = Cholesky(C);
    
    //generate the indep, unit-variance, zero-mean multivariates
    ColumnVector y(n), x(n);
    
    for (int i = 1; i <= n; i++)
        y(i) = get_gaussian(0,1);
    x = L*y;
    for (int i = 1; i <= n; i++)
        x(i) += mean(i);
    return x;
}

int main (int argc, char* argv[])
{
    ColumnVector y;
    Matrix count(100,100);
    int no_of_trials, dimension;
    
    // 2D case
    dimension = 2;
    
    sscanf (argv[1], "%d", &no_of_trials);
    ofstream pdf_data(argv[2]);
    ofstream pdf_theory(argv[3]);
    
    // The covariance matrix
    SymmetricMatrix C(2);
    C(1,1) = 1.0;
    C(1,2) = 0.5;
    C(2,1) = 0.5;
    C(2,2) = 1.0;
    
    // The mean vector
    ColumnVector mean(2);
    mean(1) = 1.0;
    mean(2) = 2.0;
    
    for (int i = 1; i <= 100; i++)
        for (int j = 1; j <= 100; j++)
            count(i,j) = 0.0;
    
    for (int i = 0; i < no_of_trials; i++) {
        y = multivariate_gaussian(dimension, C, mean);
        for (int j = 1; j <= 100; j++) {
            for (int k = 1; k <= 100; k++) {
                if ( (y(1) >= ((float) (j-52)/10)) && (y(1) < ((float) (j-51)/10)) &&
                    (y(2) >= ((float) (k-52)/10)) && (y(2) < ((float) (k-51)/10)) )
                    count(j,k)++;
            }
        }
    }
    
    for (int j = 1; j <= 100; j++) {
        for (int k = 1; k <= 100; k++) {
            if (k < 100)
                pdf_data << count(j,k)/((float) no_of_trials) << ", ";
            if (k == 100)
                pdf_data << count(j,k)/((float) no_of_trials) << endl;
        }
    }
    
    float x1, x2;
    for (int j = 1; j <= 100; j++) {
        x1 = ((float) (j-51)/10);
        for (int k = 1; k <= 100; k++) {
            x2 = ((float) (k-51)/10);
            if (k < 100)
                pdf_theory << 0.1838*exp(-(((4*x1)/3 - (2*x2)/3)*(x1 - 1) - (x2 - 2)*((2*x1)/3 - (4*x2)/3 + 2))/2)*0.01 << ", ";
            if (k == 100)
                pdf_theory << 0.1838*exp(-(((4*x1)/3 - (2*x2)/3)*(x1 - 1) - (x2 - 2)*((2*x1)/3 - (4*x2)/3 + 2))/2)*0.01 << endl;
        }
    }
}

