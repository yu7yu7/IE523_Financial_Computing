#include <iostream>
#include "binomial.h"

int main (int argc, char * const argv[]) 
{
	long n, r;
	std::cout << "What is the value of n? ";
	std::cin >> n;
	std::cout << "What is the value of r? ";
	std::cin >> r;
	std::cout << "Regular Binomial Coefficient Computation" << std::endl;
	std::cout << "(nCr) = " << Regular_Binomial(n,r) << std::endl;
    std::cout << "Recursive Binomial Coefficient Computation" << std::endl;
	std::cout << "(nCr) = " << Recursive_Binomial(n,r) << std::endl;
    return 0;
}
