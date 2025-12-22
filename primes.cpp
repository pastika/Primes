#include <cmath>
#include <ctime>
#include <iostream>
#include <cstdint>

const int NEVTS = 1000000000;
const int A = NEVTS/2;
const int SA = (int)sqrt((double)NEVTS)/2;

int main()
{
	clock_t t0, t1;
	unsigned int k, i, numPrimes = 0;
	static bool P[A];
	
	t0 = clock();	
	for(i = 0; i < A; i++) P[i] = true;
	for(k = 1; k <= SA; k++) if(P[k]) for(i = 2*k*(k+1); i < A; i += (2*k + 1)) P[i] = false;
	
	for(i = 0; i < A; i++) if(P[i]) numPrimes++;
	t1 = clock();
	
	std::cout << "I have found " << numPrimes + 1 << " primes in " << (double)(t1 - t0)/CLOCKS_PER_SEC << " seconds.\n";
}
