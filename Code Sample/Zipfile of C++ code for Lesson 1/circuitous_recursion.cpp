// Illustrative code about passing by reference in C++
#include <iostream>
using namespace std;

int main() 
{
	int c;
	if ( (c = cin.get()) != static_cast<int> ('\n')) {
		main();
		cout << static_cast<char> (c);
	}
}


