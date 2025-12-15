#include <cmath>
#include <ctime>
#include <iostream>
#include <cstdint>

const int NEVTS = 1000000000;
const int A = NEVTS/2;
const int SA = (int)sqrt((double)NEVTS)/2;
const int K = 5;

int main()
{
	clock_t t0, t1;
	unsigned int k, i, j, numPrimes = 0, M = 1;
	//static bool P[A];
	int *W, *X, *X2;
	static int S[K-1];

	t0 = clock();

	//find first k primes
	j = 0;
	for(k = 3; j < K - 1; k += 2)
	{
		for(i = 2; k%i; i++)
		{
			if(i*i >= k)
			{
				S[j++] = k;
				break;
			}
		}
	}

	for(k = 0; k < K - 1; ++k) M *= S[k];

	//calculate W
	W = new int[M];
	for(k = 0; k < M; ++k) W[k] = 1;

	for(k = 0; k < K - 1; ++k) for(i = (S[k] - 1)/2; i < M; i +=S[k]) W[i] = 0;
	
	j = M;
	for(k = M - 1; k < M; --k)
	{
		if(W[k]) 
		{
			W[k] = j - k;
			j = k;
		}
	}

	//recreate W (call in X)
	int nW = 0;
	for(i = 0; i < M; ++i) if(W[i]) nW++;
	X = new int[nW];
	j = 0;
	for(i = 0; i < nW; i++) 
	{
		X[i] = W[j];
		j += W[j];
	}
	//for(i = 0; i < nW; i++) 
	//{
	//	X2[i] = 0;
	//	for(j = 0; j < i + 1; j++) X2[i] += X[j];
	//}

	//make new P
	unsigned int tL = A*nW/M;
	for(int fgh = (tL * M) / nW; fgh < A; fgh += X[tL++%nW]);
	const unsigned int L = tL - 1;
	bool *NP = new bool[L];

	//prepair sieve
	for(k = 0; k < L; ++k) NP[k] = true;
	NP[0] = false;

	//run sieve
	int m1 = X[0], m2 = 0, m3 = 0;
	for(k = 1; m1 < SA; ++k)
	{
		if(NP[k])
		{
			m2 = m1;
			for(i = k; i < k + nW; i++)
			{
				m3 = ((2*m1+1)*(2*m2+1) - 1)/2;
				j = (m3 / M) * nW;
				for(m3 = M*(m3 / M); m3 < 2*m1*m2 + m1 + m2; ++j) m3 += X[j%nW];
				for(; j < L; j += nW*(2*m1+1)) NP[j] = false;
				m2 += X[i%nW];
			}
		}
		m1 += X[k%nW];
	}
	
	for(i = 0; i < L; i++) if(NP[i]) numPrimes++;
	t1 = clock();
		
	std::cout << "I have found " << numPrimes + 1 + K << " primes in " << (double)(t1 - t0)/CLOCKS_PER_SEC << " seconds.\n";
}
