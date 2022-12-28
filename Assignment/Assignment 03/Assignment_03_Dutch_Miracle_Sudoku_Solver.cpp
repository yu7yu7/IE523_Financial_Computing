#include <iostream>
#include "Assignment_03_Dutch_Miracle_Sudoku_Solver.h"

int main (int argc, char * const argv[]) 
{
	Sudoku x;
    x.read_puzzle(argc, argv);
    x.Solve(0,0);
	
    return 0;
}
