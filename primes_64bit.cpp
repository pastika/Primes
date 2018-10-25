#include <cmath>
#include <ctime>
#include <iostream>

const unsigned __int64 NEVTS = 10000000000ul;
unsigned __int64 A = NEVTS/2;
const unsigned __int64 SA = (unsigned __int64)sqrt((double)NEVTS)/2;

int main()
{
	clock_t t0, t1;
	unsigned __int64 k, i, numPrimes = 0;
	bool *P = new bool[A];
	
	t0 = clock();
	std::cout << "HELLO" << std::endl;
	for(i = 0; i < A; i++) P[i] = true;
	std::cout << "THERE" << std::endl;
	for(k = 1; k <= SA; k++) if(P[k]) for(i = 2*k*(k+1); i < A; i += (2*k + 1)) P[i] = false;
	t1 = clock();
	
	for(i = 0; i < A; i++) if(P[i]) numPrimes++;
	
	std::cout << "I have found " << numPrimes + 1 << " primes in " << (double)(t1 - t0)/CLOCKS_PER_SEC << " seconds.\n";
}