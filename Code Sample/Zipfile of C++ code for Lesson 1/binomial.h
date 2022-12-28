/*
 *  binomial.h
 *  Binomial Coefficient
 *
 *  Created by Ramavarapu Sreenivas on 8/31/12.
 */

#ifndef binomial
#define binomial

template <typename T>
T Recursive_Binomial(T n, T i)
{
	if (n < 0 || i < 0) 
		return 0;
	if (i > n) 
		return 0;
	if (0 == i) 
		return 1;
	if (1 == i) 
		return n;
	if (n == i) 
		return 1;
	if (n == i+1) 
		return n;
	return (Recursive_Binomial(n-1,i-1) + Recursive_Binomial(n-1,i));
}

template <typename T>
T factorial (T n)
{
	if (n < 0)
		return 0;
	if ((0 == n) || (1 == n))
		return 1;
	return (n * factorial(n-1));
}

template <typename T>
T Regular_Binomial(T n, T i)
{
	if (n < 0 || i < 0) 
		return 0;
	if (i > n) 
		return 0;
	if (0 == i) 
		return 1;
	if (1 == i) 
		return n;
	if (n == i) 
		return 1;
	if (n == i+1) 
		return n;
	return (factorial(n)/(factorial(i)*factorial(n-i)));
}
#endif