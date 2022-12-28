// Illustrative code about passing by reference in C++
#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

double my_cos(double), my_sin(double);

int main() 
{
	cout << "Computing Sine & Cosine functions recursively" << endl;
	cout << "Rads.\tmy_sin() error\tmy_cos() error" << endl;
	cout << setprecision(6);
	for (double angle = 0.1; angle <= 1.0; angle = angle + 0.1) 
		cout << angle << "\t" << my_sin(angle)-sin(angle) << "\t" << my_cos(angle)-cos(angle)<< endl;	
}

double my_sin(double theta)
{
	if ( (-0.005 < theta) && (theta < 0.005) )
		return (theta - (theta*theta*theta)/6);
	else
		return (2*my_sin(theta/2)*my_cos(theta/2));
}

double my_cos(double theta)
{
	if ( (-0.005 < theta) && (theta < 0.005) )
		return (1 - (theta*theta)/2);
	else
		return (1 - 2*my_sin(theta/2)*my_sin(theta/2));
}
