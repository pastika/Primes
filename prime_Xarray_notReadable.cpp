#include <ctime>
#include <cmath>
#include <iostream>

//Type of variable to bit pack
typedef unsigned __int16 T;
//Number of numbers to check (must be multiple of sizeof(T)*16)
const unsigned int NEVTS = 1000000000;

template< unsigned int N> class LOG2
{
public:
	static const unsigned int RESULT = 1 + LOG2<N/2>::RESULT;
};

template<> class LOG2< 2 >
{
public:
	static const unsigned int RESULT = 1;
};

//These values are automatically calculated
const unsigned int N = sizeof(T)*8;
const unsigned int A = NEVTS/(N*2);
const unsigned int SNA = (unsigned int)sqrt((double)NEVTS)/2;
const unsigned int L2N = LOG2<N>::RESULT;

int main()
{
	static T P[A];

	clock_t t0, t1;
	unsigned int numPrimes = 0, k, i;
	
	t0 = clock();
	for(i = 0; i < A; i++) P[i] = ~(T)(0);
	for(k = 1; k <= SNA; k++) if(P[k >> L2N] & (T(1) << (k & (N - 1)))) for(i = 3*k + 1; i < NEVTS/2; i += (2*k + 1)) /*if(P[i >> L2N] & (T(1) << (i & (N - 1))))*/ P[i >> L2N] &= ~(T(1) << (i & (N - 1)));
	t1 = clock();

	for(i = 0; i < NEVTS/2; i++) if(P[i >> L2N] & (T(1) << (i & (N - 1)))) numPrimes++;
	std::cout << "I have found " << numPrimes + 1 << " primes in " << (double)(t1 - t0)/CLOCKS_PER_SEC << " seconds.\n";
}

