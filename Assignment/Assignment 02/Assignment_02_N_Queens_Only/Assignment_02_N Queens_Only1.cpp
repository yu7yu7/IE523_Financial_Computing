//
// Created by LYC on 2022/09/07.
//

#include <iostream>
#include "Assignment_02_N Queens_Only1.h"

int main (int argc, char * const argv[]) 
{
    Board x;
    
    int board_size;
    sscanf (argv[1], "%d", &board_size);
	
	x.nQueens(board_size);
	
    return 0;
}
