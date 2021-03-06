#include <cmath>
#include <iostream>
#include <vector>

namespace pnlg
{
	//Type of variable to bit pack
	typedef unsigned __int64 T;
	//Number of numbers to check (must be multiple of sizeof(T)*16)
	const unsigned int NEVTS = 1000000;
	//Number of arrays of size NEVTS to search through (must be less than NEVTS)
	const unsigned int NN = 1000;

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

	inline unsigned int mymod(unsigned int a, unsigned int b, unsigned int c);

	std::vector<unsigned int> getPrimeList(unsigned int num2check)
	{
		unsigned int n, k, i, j, svmk, numPrimes = 0, l;
		static T P[A], Q[A];
		T im;
		std::vector<unsigned int> primeList;
		primeList.reserve((unsigned int)(1.2*num2check/log((double)num2check)));
		for(i = 0; i < A; i++) Q[i] = ~(T)(0);
		for(k = 1; k <= SNA; k++)
		{
			if(Q[k >> L2N] & (T(1) << (k & (N - 1)))) 
			{
				if(k < N/2)
				{
					l = (j = 2*k+1) - N%j;
					svmk = k - l;
					im = T(1);
					for(i = 0; i < N; i += j) im |= im << j;
					for(i = 0; i < A; i++) Q[i] &= ~(im << (((svmk += l) < j)?svmk:(svmk -= j)));
					Q[k >> L2N] |= T(1) << k;
				}
				else for(i = 2*k*(k+1); i < NEVTS/2; i += (2*k + 1)) Q[i >> L2N] &= ~(T(1) << (i & (N - 1)));
			}
		}
	
		for(i = 0; i < NEVTS/2; i++) 
			if(Q[i >> L2N] & (T(1) << (i & (N - 1)))) 
			{
				numPrimes++;
				primeList.push_back(2*i + 1);
				if(2*i + 1 > num2check) return primeList;
			}
	
		for(n = 1; n < NN; n++)
		{
			for(i = 0; i < A; i++) P[i] = ~(T)(0);
			//for(k = 2*(j = 1) + 1; k <= NEVTS/2 && (float)k/(n+1) < NEVTS/(float)k; k = 2*(++j) + 1) if(Q[j >> L2N] & (T(1) << (j & (N - 1)))) for(i = ((svmk = mymod(NEVTS, n, k)) & 1)?(k - (svmk+1)/2):((k - svmk-1)/2); i < NEVTS/2; i += k) P[i >> L2N] &= ~(T(1) << (i & (N - 1)));
			for(k = 2*(j = 1) + 1; k <= NEVTS/2 && (float)k/(n+1) < NEVTS/(float)k; k = 2*(++j) + 1) 
			{
				if(Q[j >> L2N] & (T(1) << (j & (N - 1))))
				{
					if(k < N)
					{
						l = k - N%k;
						svmk = (((svmk = mymod(NEVTS, n, k)) & 1)?(k - (svmk+1)/2):((k - svmk-1)/2)) - l;
						im = T(1);
						for(i = 0; i < N; i += k) im |= im << k;
						for(i = 0; i < A; i++) P[i] &= ~(im << (((svmk += l) < k)?svmk:(svmk -= k)));
					}
					else for(i = ((svmk = mymod(NEVTS, n, k)) & 1)?(k - (svmk+1)/2):((k - svmk-1)/2); i < NEVTS/2; i += k) /*if(P[i >> L2N] & (T(1) << (i & (N - 1))))*/ P[i >> L2N] &= ~(T(1) << (i & (N - 1)));
				}
			}
			for(i = 0; i < NEVTS/2; i++) 
				if(P[i >> L2N] & (T(1) << (i & (N - 1)))) 
				{
					numPrimes++;
					primeList.push_back(n*NEVTS + 2*i + 1);
					if(n*NEVTS + 2*i + 1 > num2check) return primeList;
				}
		}

		return primeList;
	}

	unsigned int mymod(unsigned int a, unsigned int b, unsigned int c)
	{
		unsigned int d = a%c, f = b%c, r = 0, rlast = 0, inc = d?(c/d):1;
		for(unsigned int i = 0; i < f; (f - i < inc)?(i++):(i += inc)) 
		{
			if(((f - i < inc)?(r += d):(r += inc*d)) > c || rlast > r) r -= c;
			rlast = r;
		}
		return r;
	}

}