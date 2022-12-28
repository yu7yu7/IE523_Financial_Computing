// Recursion in C++ illustration -- Tower of Hanoi Problem
// Written by Prof. Sreenivas for IE523: Financial Computing
#include <iostream>
using namespace std;

// recursive function that computes the required moves
void move(int n, char source, char destination, char intermediate) 
{ 
	if (n > 0) { 
		move (n-1, source, intermediate, destination); 
		cout << "Move the top disk from peg " << source << " to peg " << destination << endl;
		move (n-1, intermediate, destination, source); 
	} 
	/* if n == 0, do nothing */ 
	return;
} 

int main () 
{ 
	int no_of_disks; 
	cout << "Enter the #disks in peg A: ";
	cin >> no_of_disks;
	
	if (no_of_disks <= 10) 
		move(no_of_disks, 'A', 'C', 'B'); 
	else { 
		cout << "This is going to take too long... exiting" << endl;
		exit (0); 
	} 
} 




