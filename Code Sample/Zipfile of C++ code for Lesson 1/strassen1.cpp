// (Tons of) Recursion in C++ illustration -- Strassen's Method of 
// Matrix Multiplication 
// Written by Prof. Sreenivas for IE523: Financial Computing
// See section 3.2 of Lesson 2 of my notes for the relevant details.

#include <iostream>
#include <cmath>
using namespace std;

// This is the standard matrix multiplication procedure
// It takes as input two n x n matrices called matrix_1 and matrix_2 
// computes the product and puts the product in an n x n matrix called result.
// Note that the size n is passed as a parameter to the routine, which
// is then used to reserve appropriate amount of memory (contrast this with 
// the process where a very very large storage space is reserved for to 
// store the matrix and you end up using just a small portion of it in 
// the specific computation). 
double **regular_multiply(double **matrix_1, double **matrix_2, int n)
{
	int i, j, k; 
	
	// creating the necessary storage for the product matrix
	// note, the size of this matrix is dynamic (i.e. it is decided by n, 
	// which is an input parameter.
	double **result = new double*[n];
	for (i = 0; i < n; i++)
		result[i] = new double[n];
	
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			result[i][j] = 0.0;
			for (k = 0; k < n; k++) 
				result[i][j] = result[i][j] + matrix_1[i][k]*matrix_2[k][j];
		}
	}
	// return the pointer to the product matrix
	return (result);
}

// Routine to add two n x n matrices and return the result
// Note the size n is a variable (see the blurb in the comments that precede
// the "regular_multiply" routine listed above.
double **matrix_add(double **m1, double **m2, int n)

{
	int i, j;
	// creating the necessary storage for the product matrix
	// note, the size of this matrix is dynamic (i.e. it is decided by n, 
	// which is an input parameter.
	double **result = new double*[n];
	for (i = 0; i < n; i++)
		result[i] = new double[n];
	
	for (i = 0; i < n; i++) 
		for (j = 0; j < n; j++) 
			result[i][j] = m1[i][j] + m2[i][j];
	
	// return the pointer to the product matrix
	return (result);
}

// Routine to subtract two n x n matrices and return the result
// Note the size n is a variable (see the blurb in the comments that precede
// the "regular_multiply" routine listed above.
double **matrix_subtract(double **m1, double **m2, int n)

{
	int i, j;
	// creating the necessary storage for the product matrix
	// note, the size of this matrix is dynamic (i.e. it is decided by n, 
	// which is an input parameter.
	double **result = new double*[n];
	for (i = 0; i < n; i++)
		result[i] = new double[n];
	
	for (i = 0; i < n; i++) 
		for (j = 0; j < n; j++) 
			result[i][j] = m1[i][j] - m2[i][j];
	
	// return the pointer to the product matrix
	return (result);
}

// This recursive procedure implements the Strassen Multiplication Routine
// It takes as input two (n x n) matrices and computes their product using
// a recursive procedure and puts the product in a matrix called result
// It is important to note that the procedure supposes that n = 2^k for 
// some k (i.e. you can keep dividing n by 2; you will keep getting an even
// number)
double **strassen_multiply(double **matrix_1, double **matrix_2, int n)

{
	// some integer variables that I need for the for-loops and temporary work
	int i, j, k;
	
	// creating the necessary storage for the product matrix
	// note, the size of this matrix is dynamic (i.e. it is decided by n, 
	// which is an input parameter.
	j = (int) pow(2.0, (double) n);
	double **result = new double*[j];
	for (i = 0; i < j; i++)
		result[i] = new double[j];
	
	// Let us pull some teeth now!!!
	
	if (n == 1) {
		// this is my exit condition for the recursion; you can exit with a 
		// regular multiply at some other value of n (say n = 2) or something
		// you will usually figure that out after calibrating your implementation.
		k = (int) pow(2.0, (double) n);
		result = regular_multiply(matrix_1, matrix_2, k);
	}
	else {
		
		k = (int) pow(2.0, (double) n-1);
		// creating the appropriate storage space 
		
		// the 4 (n/2 x n/2) submatrices of the two matrices to be multiplied 
		// are identified as **new_m1_11, **new_m1_12, **new_m1_21, **new_m1_22;
		// **new_m2_11, **new_m2_12, **new_m2_21, **new_m2_22;
		double **new_m1_11 = new double*[k];
		double **new_m1_12 = new double*[k];
		double **new_m1_21 = new double*[k];
		double **new_m1_22 = new double*[k];

		double **new_m2_11 = new double*[k];
		double **new_m2_12 = new double*[k];
		double **new_m2_21 = new double*[k];
		double **new_m2_22 = new double*[k];

		// I need these 4 submatrices for the final product 
		// **new_r_11, **new_r_12, **new_r_21, **new_r_22;
		double **new_r_11 = new double*[k];
		double **new_r_12 = new double*[k];
		double **new_r_21 = new double*[k];
		double **new_r_22 = new double*[k];

		// storage variables to make my life easy for intermediate steps
		// **temp1, **temp2;
		double **temp1 = new double*[k];
		double **temp2 = new double*[k];
		
		// See my notes to figure out what these seven alpha matrices are
		// **alpha1, **alpha2, **alpha3, **alpha4, **alpha5, **alpha6, **alpha7;
		double **alpha1 = new double*[k];
		double **alpha2 = new double*[k];
		double **alpha3 = new double*[k];
		double **alpha4 = new double*[k];
		double **alpha5 = new double*[k];
		double **alpha6 = new double*[k];
		double **alpha7 = new double*[k];
		
		for (i = 0; i < k; i++) {
			
			new_m1_11[i] = new double[k];
			new_m1_12[i] = new double[k];
			new_m1_21[i] = new double[k];
			new_m1_22[i] = new double[k];

			new_m2_11[i] = new double[k];
			new_m2_12[i] = new double[k];
			new_m2_21[i] = new double[k];
			new_m2_22[i] = new double[k];
			
			new_r_11[i] = new double[k];
			new_r_12[i] = new double[k];
			new_r_21[i] = new double[k];
			new_r_22[i] = new double[k];
			
			temp1[i] = new double[k];
			temp2[i] = new double[k];
			
			alpha1[i] = new double[k];
			alpha2[i] = new double[k];
			alpha3[i] = new double[k];
			alpha4[i] = new double[k];
			alpha5[i] = new double[k];
			alpha6[i] = new double[k];
			alpha7[i] = new double[k];

		}
		
		// filling the values of the newly created matrices
		for (i = 0; i < k; i++) {
			for (j = 0; j < k; j++) {
				
				new_m1_11[i][j] = matrix_1[i][j];
				new_m1_12[i][j] = matrix_1[i][j+k];
				new_m1_21[i][j] = matrix_1[i+k][j];
				new_m1_22[i][j] = matrix_1[i+k][j+k];
				
				new_m2_11[i][j] = matrix_2[i][j];
				new_m2_12[i][j] = matrix_2[i][j+k];
				new_m2_21[i][j] = matrix_2[i+k][j];
				new_m2_22[i][j] = matrix_2[i+k][j+k];
			}
		}
		
		temp1 = matrix_add(new_m1_11, new_m1_22, k);
		temp2 = matrix_add(new_m2_11, new_m2_22, k);
		alpha1 = strassen_multiply(temp1, temp2, n-1);
		
		temp1 = matrix_subtract(new_m2_21, new_m2_11, k);
		alpha2 = strassen_multiply(new_m1_22, temp1, n-1);
		
		temp1 = matrix_subtract(new_m1_12, new_m1_22, k);
		temp2 = matrix_add(new_m2_21, new_m2_22, k);
		alpha3 = strassen_multiply(temp1, temp2, n-1);
		
		temp1 = matrix_add(new_m1_11, new_m1_12, k);
		alpha4 = strassen_multiply(temp1, new_m2_22, n-1);
		
		temp1 = matrix_add(new_m1_21, new_m1_22, k);
		alpha5 = strassen_multiply(temp1, new_m2_11, n-1);
		
		temp1 = matrix_subtract(new_m2_12, new_m2_22, k);
		alpha6 = strassen_multiply(new_m1_11, temp1, n-1);
		
		temp1 = matrix_subtract(new_m1_21, new_m1_11, k);
		temp2 = matrix_add(new_m2_11, new_m2_12, k);
		alpha7 = strassen_multiply(temp1, temp2, n-1);
		
		temp1 = matrix_add(alpha1, alpha2, k);
		temp2 = matrix_subtract(alpha3, alpha4, k);
		new_r_11 = matrix_add(temp1, temp2, k);
		
		new_r_12 = matrix_add(alpha4, alpha6, k);
		
		new_r_21 = matrix_add(alpha2, alpha5, k);
		
		temp1 = matrix_subtract(alpha1, alpha5, k);
		temp2 = matrix_add(alpha6, alpha7, k);
		new_r_22 = matrix_add(temp1, temp2, k);
		
		for (i = 0; i < k; i++) {
			for (j = 0; j < k; j++) {
				result[i][j] = new_r_11[i][j];
				result[i][j+k] = new_r_12[i][j];
				result[i+k][j] = new_r_21[i][j];
				result[i+k][j+k] = new_r_22[i][j];
			}
		}
	}

	// return the product matrix
	return (result);
}

// The following routine returns a uniformly distributed random number in 
// the unit-interval... we will cover more of how this works later...
double get_uniform()
{
	return (((double) random())/(pow(2.0, 31.0)-1.0));
}

// This routine fills a n x n matrix, with random entries
void create_random_matrix(double **matrix, int n)

{
	int i, j;
	
	for (i = 0; i < n; i++) 
		for (j = 0; j < n; j++) 
			matrix[i][j] = 100.0*get_uniform();
	
}	

// I am not using this routine -- printing a (2^k x 2^k) matrix for large 
// k can be purposeless... nonetheless, I thought I would leave it in here
// if you need to debug/check things.
void print_matrix(double **m, int n)
{
	int i, j;
	
	cout << endl;
	cout << "---------------";
	
	for (i = 0; i < n; i++) {
		cout << endl;;
		for (j = 0; j < n; j++) 
			cout << m[i][j]; 
	}
	cout << endl;
	cout << "---------------";	
}

// In this main routine I am generating two random (2^k x 2^k) matrices, while
// k varies from 2 to 8 (i.e. two (2 x 2), (4 x 4), (8 x 8) ... (256 x 256) randomly
// generated matrices).  I compute the product of these two (randomly generated) 
// matrices by Strassen's method and regular matrix multiplication. I present the
// cummulative errors between each member of the matrix as a measure of accurracy.

int main ()
{
	int my_index, k;  
	double error = 0.0;
	
	srand( (unsigned)time( NULL ) );

	for (my_index = 1; my_index < 9; my_index++) {
		
		// creating the space for four matrices to be used for testing
		// if Strassen's method works
		k = (int) pow(2.0, (double) my_index);
		
		double **first_matrix = new double*[k];
		double **second_matrix = new double*[k];
		double **third_matrix = new double*[k];
		double **fourth_matrix = new double*[k];
 
		for (int i = 0; i < k; i++) {
			first_matrix[i] = new double[k];
			second_matrix[i] = new double[k];
			third_matrix[i] = new double[k];
			fourth_matrix[i] = new double[k];
		}

		// fill the first and second matrix with random entries
		create_random_matrix(first_matrix, k);
		create_random_matrix(second_matrix, k);

		third_matrix = strassen_multiply(first_matrix, second_matrix, my_index);

		fourth_matrix = regular_multiply(first_matrix, second_matrix, k);

		for (int j = 0; j < k; j++) 
			for (int i = 0; i < k; i++) 
				error += fabs(third_matrix[i][j] - fourth_matrix[i][j]);
		
		cout << "Matrix Size = (" << k << " x " << k << "); Absolute Error = " << error << endl;
		
		// free the space used up by the matrices
		delete [] first_matrix, second_matrix, third_matrix, fourth_matrix;
	}
	
}



