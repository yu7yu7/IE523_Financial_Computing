// Matrix Exponentiation by Repeated Squaring
// We need NEWMAT library for this... this would mean you have
// to set all the include, library and linker flags according to
// your compiler.
// Written by Prof. Sreenivas for IE523: Financial Computing
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include "/Users/sreenivas/Documents/Courses/IE523/newmat10/newmat.h"

float get_uniform()
{
	return (((float) random())/(pow(2.0, 31.0)-1.0));
}

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

int main (int argc, char* argv[])
{
	int exponent, dimension;
	clock_t time_before, time_after;
	float diff;
	
	sscanf (argv[1], "%d", &exponent);
	sscanf (argv[2], "%d", &dimension);
	
	cout << "Comparison of brute-force matrix exponentiation vs. Repeated Squaring" << endl;
	cout << "Exponent = " << exponent << endl;
	cout << "Matrix Dimension = " << dimension << endl;

	Matrix A(dimension,dimension), B(dimension,dimension), C(dimension,dimension);

	for(int i = 1; i <= dimension; i++)
		for (int j = 1; j <= dimension; j++)
			A(i,j) = 10*(get_uniform()-0.5);
	
	time_before = clock(); // recording time before repeated squaring algorithm starts
	B = repeated_squaring(A, exponent, dimension);
	time_after = clock(); // recording time after repeated squaring algorithm finishes
	diff = ((float) time_after - (float) time_before);
	cout << "---------------" << endl;
	cout << "Repeated Squaring took " << diff/CLOCKS_PER_SEC << " seconds to complete" << endl;

	for (int i = 1; i <= dimension; i++) {
		cout << endl;
		for (int j = 1; j <= dimension; j++) 
			cout << B(i,j) << " ";
	}
	cout << endl;
	
	time_before = clock();
	C = A;
	for (int i = 1; i < exponent; i++) 
		C = A*C;
	time_after = clock();
	diff = ((float) time_after - (float) time_before);
	cout << "---------------" << endl;
	cout << "Brute Force Multiplication took " << diff/CLOCKS_PER_SEC << " seconds to complete" << endl;
	
	for (int i = 1; i <= dimension; i++) {
		cout << endl;
		for (int j = 1; j <= dimension; j++) 
			cout << C(i,j) << " ";
	}
	cout << endl;	
}
	
	





