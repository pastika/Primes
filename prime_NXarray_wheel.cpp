#include <cmath>
#include <ctime>
#include <iostream>

//Type of variable to bit pack
typedef unsigned __int64 T;
//Type of unsigned int
typedef unsigned int uint;
//Type of unsigned int 64
typedef unsigned __int64 uint64;
//Number of numbers to check (must be multiple of sizeof(T)*16)
const uint64 NEVTS = 10000000;
//Number of arrays of size NEVTS to search through (must be less than NEVTS)
const uint64 NN = 100;
//Wheel number
const uint64 K = 3;

//Calculates log base 2 (this really finds how many times 2 is a factor of N)
template<unsigned int N, unsigned int G = (N & 1)> class LOG2
{
public:
	static const unsigned int RESULT = 1 + LOG2<N/2, (N/2) & 1>::RESULT;
};

template<unsigned int N> class LOG2<N, 1>
{
public:
	static const unsigned int RESULT = 0;
};

//Calculates 2^N
template<unsigned int N> class POW2
{
public:
	static const unsigned int RESULT = 2*POW2<N - 1>::RESULT;
};

template<> class POW2<0>
{
public:
	static const unsigned int RESULT = 1;
};

//Calculates the max of N and M
template<unsigned int N, unsigned int M, bool P = (N > M)> class MAX
{
public:
	static const unsigned int RESULT = N;
};

template<unsigned int N, unsigned int M> class MAX<N, M, false>
{
public:
	static const unsigned int RESULT = M;
};

//Calculates the LCM w.r.t. base 2 (as M is always 2^X)
template<unsigned int N, unsigned int M, unsigned int P = 0> class LCM2
{
public:
	static const unsigned int RESULT = M*N*POW2<P>::RESULT;
};

template<unsigned int N, unsigned int M> class LCM2<N, M, 0>
{
public:
	static const unsigned int RESULT = LCM2<N/POW2<LOG2<N>::RESULT>::RESULT, M/POW2<LOG2<M>::RESULT>::RESULT, MAX<LOG2<N>::RESULT, LOG2<M>::RESULT >::RESULT>::RESULT;
};

//Returns N if it is prime, 0 otherwise
template<unsigned int N, unsigned int M = N/2> class ISPRIME
{
public:
	static const unsigned int RESULT = (N % M > 0) * ISPRIME<N, M - 1>::RESULT;
};

template<unsigned int N> class ISPRIME<N, 1>
{
public:
	static const unsigned int RESULT = N;
};

//Calculates the next prime 
template<unsigned int N, bool P = false> class NEXTPRIME
{
public:
	static const unsigned int RESULT = NEXTPRIME<N + 1, (ISPRIME<N + 1>::RESULT > 0) >::RESULT;
};

template<unsigned int N> class NEXTPRIME<N, true>
{
public:
	static const unsigned int RESULT = N;
};

//Returns product of first N primes excluding 2
template<unsigned int N, unsigned int M = 3> class PRIMEPROD
{
public:
	static const unsigned int RESULT = M*PRIMEPROD<N - 1, NEXTPRIME<M>::RESULT>::RESULT;
};

template<unsigned int M> class PRIMEPROD<0, M>
{
public:
	static const unsigned int RESULT = 1;
};

//Calculate the length of the wheel array X
template<unsigned int N, unsigned int P = 2> class LENX
{
public:
	static const unsigned int RESULT = (P - 1)*LENX<N - 1, NEXTPRIME<P>::RESULT>::RESULT;
};

template<unsigned int P> class LENX<0, P>
{
public:
	static const unsigned int RESULT = P - 1;
};

//These values are automatically calculated
const uint64 N = sizeof(T)*8;
const uint64 A = NEVTS/(N*2);
const uint64 SNA = (unsigned int)sqrt((double)NEVTS)/2;
const uint64 L2N = LOG2<N>::RESULT;
const uint64 M = PRIMEPROD<K - 1>::RESULT;
const uint64 nW = LENX<K - 1>::RESULT;
const uint64 L = (A/M)*nW + (nW*(A%M))/M + 1;
const uint64 KT = (uint64)sqrt(double(N));

inline int sparse_ones_bitcount(T n);

int main()
{
    clock_t t0, t1;
	uint n, k, i, j, svmk, numPrimes = 0, l;
	uint64 m1 = 0, m2 = 0, m3 = 0, m = 0;
	uint W[M], X[nW];
	static T P[L], Q[L];
	T im;
		
	t0 = clock();

	int S[K-1];
	//find first k primes
	for(k = 3, j = 0; j < K - 1; k += 2) for(i = 2; k%i && (j == 0 || S[j-1] != k); i++) if(i*i >= k) S[j++] = k;
	
	//calculate W and X
	for(k = 0; k < M; ++k) W[k] = 1;
	for(k = 0; k < K - 1; ++k) for(i = (S[k] - 1)/2; i < M; i += S[k]) W[i] = 0;
	for(k = (j = M) - 1; k < M; --k)
	{
		if(W[k]) 
		{
			W[k] = j - k;
			j = k;
		}
	}
	for(i = j = 0; i < nW; i++) j += X[i] = W[j];
	
	for(k = 0; k < L; ++k) Q[k] = ~(T)(0);
	for(k = m1 = 0; k < KT; m1 += X[k++%nW]);
	const uint M2M = (2*m1+1)*(LCM2<nW, N>::RESULT/N);
	T *B = new T[M2M];
	m1 = X[0];
	for(k = 1; k < KT; ++k)
	{
		if(Q[k >> L2N] & (T(1) << (k & (N - 1)))) 
		{
			//calculate B
			const uint M2 = (2*m1+1)*(LCM2<nW, N>::RESULT/N);
			for(i = 0; i < M2; ++i) B[i] = ~(T)(0);
			m2 = 2*m1+1;
			for(i = k; i < M2*N; ++i)
			{
				if(m2 == 2*m1+1) 
				{
					B[i >> L2N] &= ~(T(1) << (i & (N - 1)));
					m2 -= 2*m1+1;
				}
				else if(m2 > 2*m1+1) m2 -= 2*m1+1;
				m2+=2*X[i%nW];
			}
			for(i = 0; i < L; i++) Q[i] &= B[i%M2];
			Q[k >> L2N] |= T(1) << (k & (N - 1));
		}
		m1 += X[k%nW];
	}

	for(k = KT; k <= SNA*nW/M + nW; k++)
	{
		if(Q[k >> L2N] & (T(1) << (k & (N - 1)))) 
		{
			m2 = m1;
			for(i = 0; i < nW; i++)
			{
				m3 = ((2*m1+1)*(2*m2+1) - 1)/2;
				j = (m3 / M) * nW;
				for(m3 = M*(m3 / M); m3 < 2*m1*m2 + m1 + m2; ++j) m3 += X[j%nW];
				for(; j < L*N; j += nW*(2*m1+1)) Q[j >> L2N] &= ~(T(1) << (j & (N - 1)));
				m2 += X[(k+i)%nW];
			}
		}
		m1 += X[k%nW];
	}

	for(i = 0; i < L; i++) numPrimes += sparse_ones_bitcount(Q[i]);
	//t1 = clock();
	//std::cout << "I have found " << numPrimes + K << " primes in " << (double)(t1 - t0)/CLOCKS_PER_SEC << " seconds.\n";
	
	for(n = 1; n < NN; n++)
	{
		for(i = 0; i < L; i++) P[i] = ~(T)(0);
		m1 = X[0];
		for(k = 1; k < KT; ++k)
		{
			if(Q[k >> L2N] & (T(1) << (k & (N - 1)))) 
			{
				//calculate B
				const uint M2 = (2*m1+1)*(LCM2<nW, N>::RESULT/N);
				for(i = 0; i < M2; ++i) B[i] = ~(T)(0);
				l = (nW*(2*m1+1)) - (n*N*L)%(nW*(2*m1+1));
				m2 = 2*m1+1;
				for(i = k; i < M2*N; ++i)
				{
					if(m2 == 2*m1+1)
					{
						B[(i + l)%(N*M2) >> L2N] &= ~(T(1) << ((i + l)%(N*M2) & (N - 1)));
						m2 -= 2*m1+1;
					}
					else if(m2 > 2*m1+1) m2 -= 2*m1+1;
					m2 += 2*X[i%nW];
				}
				for(i = 0; i < L; i++) P[i] &= B[i%M2];
			}
			m1 += X[k%nW];
		}
		
		for(k = KT; (2*m1+1) < sqrt(double(n+1)*((double(2)*double(N)*double(L)*double(M))/double(nW))); k++)
		{
			if(Q[k >> L2N] & (T(1) << (k & (N - 1)))) 
			{
				m2 = m1;
				l = (nW*(2*m1+1)) - (uint64(n)*N*L)%(nW*(2*m1+1));
				for(i = 0; i < nW; i++)
				{
					m3 = ((2*m1+1)*(2*m2+1) - 1)/2;
					m = (m3 / M) * nW;
					for(m3 = M*(m3 / M); m3 < 2*m1*m2 + m1 + m2; ++m) m3 += X[m%nW];
					if(m < n*N*L)           for(j = (m+l)%(nW*(2*m1+1)); j < L*N; j += nW*(2*m1+1)) P[j >> L2N] &= ~(T(1) << (j & (N - 1)));
					else if (m < (n+1)*N*L) for(j = m%(N*L)            ; j < L*N; j += nW*(2*m1+1)) P[j >> L2N] &= ~(T(1) << (j & (N - 1)));
					m2 += X[(k+i)%nW];
				}
			}

			m1 += X[k%nW];
		}
		for(i = 0; i < L; i++) numPrimes += sparse_ones_bitcount(P[i]);

		//t1 = clock();
	    //std::cout << "I have found " << numPrimes + K << " primes in " << (double)(t1 - t0)/CLOCKS_PER_SEC << " seconds.\n";
	}
	t1 = clock();
	std::cout << "I have found " << numPrimes + K << " primes in " << (double)(t1 - t0)/CLOCKS_PER_SEC << " seconds.\n";

	delete [] B;

	for(m1 = M*((NN*N*L*M/nW)/M), m = nW*((NN*N*L*M/nW)/M); m < ((NN*N*L)); ++m) m1 += X[m%nW];
	for(m = N*L - 1; 2*m1+1 > NEVTS*NN; m1 -= X[((m--)+(NN-1)*N*L)%nW]) if(((NN == 1)?Q:P)[(m%(N*L)) >> L2N] & (T(1) << ((m%(N*L)) & (N - 1)))) numPrimes--;

	t1 = clock();
	std::cout << "I have found " << numPrimes + K << " primes in " << (double)(t1 - t0)/CLOCKS_PER_SEC << " seconds.\n";
}

int sparse_ones_bitcount(T n)
{
    int count=0 ;
    while (n)
    {
        count++ ;
        n &= (n - 1) ;     
    }
    return count ;
}

