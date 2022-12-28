// Generic Markov Chain Solver
// Written by Prof. Sreenivas
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include "newmat.h"
#include "newmatio.h"
#include "newmatap.h" 
using namespace std;

#define ERROR 1.0e-6

double get_uniform()
{
	return (((double) random())/(pow(2.0, 31.0)-1.0));
}

// The matrix A should be a row-stochastic matrix.  That is, 
// entries in each row should sum to 1.  This function computes
// the left-eigenvector of A that corresponds to the eigenvalue 
// of 1 by an iterative (i.e. repeated multiplication of A) process.
RowVector compute_left_eigenvector(Matrix A)
{
	int n = A.nrows();
	double error = 0;
	
	Matrix A_n = A;
	
	for (int i = 2; i <= n; i++)
		for (int j = 1; j <= n; j++)
			error += abs(A_n(i,j) - A_n(1,j));
	
	while (abs(error) > ERROR)
	{
		A_n = A_n * A;
		
		error = 0;
		for (int i = 2; i <= n; i++)
			for (int j = 1; j <= n; j++)
				error += abs(A_n(i,j) - A_n(1,j));
	}
	
	return (A_n.row(1));
}

// This function figures out what the next state should be for a given 
// uniformly distributed rv x and a probability matrix x
int what_is_the_next_state(Matrix A, int current_state, double x)
{
	double divisions_of_the_unit_interval[A.nrows()+1];
	int next_state;
	
	divisions_of_the_unit_interval[0] = 0.0;
	divisions_of_the_unit_interval[A.nrows()] = 1.0;
	for (int i = 1; i < A.nrows(); i++) 
		divisions_of_the_unit_interval[i] = divisions_of_the_unit_interval[i-1] + 
		A(current_state,i);
	
	for (int i = 1; i <= A.nrows(); i++)
		if ((x > divisions_of_the_unit_interval[i-1]) &&
			(x <= divisions_of_the_unit_interval[i]))
			next_state = i;
	
	return (next_state);
}

// This function simulates the process of jumping around the states of a discrete
// time Markov Chain for "no_of_iterations" and then presents the visit ratio of each
// state as the output, at the end of the simulation

RowVector verification_by_simulation(Matrix A, int no_of_iterations)
{
	// pick a random state to start the simulation
	int current_state = (((double) A.nrows()) * get_uniform()) + 1;
	double x;
	
	RowVector no_of_visits(A.nrows());
	
	// Initializing the vector that counts the number of visits to each state
	for (int i = 1; i <= A.nrows(); i++) 
	{
		no_of_visits(i) = 0;
	}
	
	// simulating no_of_iterations many jumps/visits
	for (int i = 0; i < no_of_iterations; i++) 
	{
		no_of_visits(current_state)++;
		x = get_uniform();
		current_state = what_is_the_next_state(A, current_state, x);
	}

	RowVector visit_ratio(A.nrows());
	for (int i = 1; i <= A.nrows(); i++)
		visit_ratio(i) = ((double) no_of_visits(i))/((double) no_of_iterations);
	
	return (visit_ratio);
	
}

int main (int argc, char* argv[])
{
	int size, no_of_iterations;
	sscanf(argv[1], "%d", &size);
	sscanf(argv[2], "%d", &no_of_iterations);
	ifstream input_file(argv[3]);
	
	Matrix A(size, size);
	
	if (input_file.is_open()) 
	{
		for (int i = 1; i <= size; i++) 
		{
			for (int j = 1; j <= size; j++) 
			{
				double value_just_read;
				input_file >> value_just_read;
				A(i,j) = value_just_read;
			}
		}
		
		RowVector eigenvector = compute_left_eigenvector(A);
		
		cout << "----------------------------------------" << endl;
		cout << "Generic Markov Chain Solver" << endl << endl;
		cout << "Probability Matrix (P): " << endl;
		cout << A << endl;
		cout << "Steady State Probability Vector (Left Eigenvector x): " << endl;
		cout << eigenvector << endl;
		cout << "Verification (of Left Eigenvector; ie. x*A): " << endl;
		cout << eigenvector*A << endl;
		cout << "Visit Ratio for each state after " << no_of_iterations << " many moves" << endl;
		eigenvector = verification_by_simulation(A, no_of_iterations);
		cout << eigenvector << endl;
		cout << "----------------------------------------" << endl;
	}
	else 
	{
		// I need all this stuff to print out the name of the PWD
		char *path=NULL;
		size_t size;
		path = getcwd(path, size);
		cout << "Input file: " << argv[3] << " does not exist in " << path << endl;		
	}
	
}