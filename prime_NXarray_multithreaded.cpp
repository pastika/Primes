#include <cmath>
#include <ctime>
#include <iostream>
#include <thread>
#include <cstdint>
#include <chrono>

//Type of variable to bit pack
typedef std::uint32_t T;
//Number of numbers to check (must be multiple of sizeof(T)*16)
const unsigned int NEVTS = 100000000;
//Number of arrays of size NEVTS to search through
const unsigned int NN = 100;
//number of threads to use
const unsigned int NTHREADS = 10;

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

inline unsigned int mymod(const unsigned int a, const unsigned int b, const unsigned int c);
void calcPrimes_thread(const unsigned int nstart, unsigned int* nprimes);

static T Q[A];

int main()
{
    std::chrono::steady_clock::time_point t0, t1;
    unsigned int k, i, numPrimes = 0, nppt[NTHREADS];
    std::thread *pts[NTHREADS];
	
    t0 = std::chrono::steady_clock::now();
    for(i = 0; i < A; i++) Q[i] = ~(T)(0);
    for(k = 1; k <= SNA; k++) if(Q[k >> L2N] & (T(1) << (k & (N - 1)))) for(i = 2*k*(k+1); i < NEVTS/2; i += (2*k + 1)) Q[i >> L2N] &= ~(T(1) << (i & (N - 1)));
    for(i = 0; i < NEVTS/2; i++) if(Q[i >> L2N] & (T(1) << (i & (N - 1)))) numPrimes++;

    for(unsigned int i = 0; i < NTHREADS; i++) pts[i] = new std::thread(calcPrimes_thread, i + 1, nppt + i); 
    for(unsigned int i = 0; i < NTHREADS; i++) 
    {
        pts[i]->join();
        numPrimes += nppt[i];
    }

    t1 = std::chrono::steady_clock::now();
    std::chrono::duration<double, std::milli> elapsed = t1 - t0;
	
    std::cout << "I have found " << numPrimes + 1 << " primes in " << static_cast<double>(elapsed.count())/1000 << " seconds.\n";
}

unsigned int mymod(const unsigned int a, const unsigned int b, const unsigned int c)
{
    const unsigned int d = a%c, f = b%c, inc = d?(c/d):1;
    unsigned int rlast = 0, r = 0;
    for(unsigned int i = 0; i < f; (f - i < inc)?(i++):(i += inc)) 
    {
        if(((f - i < inc)?(r += d):(r += inc*d)) > c || rlast > r) r -= c;
        rlast = r;
    }
    return r;
}

void calcPrimes_thread(const unsigned int nstart, unsigned int* nprimes)
{
    std::chrono::steady_clock::time_point t0, t1;

    T *P = new T[A];
    unsigned int n, k, i, j, svmk, im, l;
    *nprimes = 0;
    t0 = std::chrono::steady_clock::now();
    for(n = nstart; n < NN; n += NTHREADS)
    {
        for(i = 0; i < A; i++) P[i] = ~(T)(0);
        for(k = 2*(j = 1) + 1; k <= NEVTS/2 && (float)k/(n+1) < NEVTS/(float)k; k = 2*(++j) + 1) 
        {
            if(Q[j >> L2N] & (T(1) << (j & (N - 1))))
            {
                if(k < N)
                {
                    l = k - N%k;
                    svmk = (((svmk = mymod(NEVTS, n, k)) & 1)?(k - (svmk+1)/2):((k - svmk-1)/2)) - l;
                    im = 1;
                    for(i = 0; i < N; i += k) im |= im << k;
                    for(i = 0; i < A; i++) P[i] &= ~(im << (((svmk += l) < k)?svmk:(svmk -= k)));
                }
                else for(i = ((svmk = mymod(NEVTS, n, k)) & 1)?(k - (svmk+1)/2):((k - svmk-1)/2); i < NEVTS/2; i += k) if(P[i >> L2N] & (T(1) << (i & (N - 1)))) P[i >> L2N] &= ~(T(1) << (i & (N - 1)));
            }
        }
        for(i = 0; i < NEVTS/2; i++) if(P[i >> L2N] & (T(1) << (i & (N - 1)))) (*nprimes)++;
    }
    t1 = std::chrono::steady_clock::now();
    std::chrono::duration<double, std::milli> elapsed = t1 - t0;
    std::cout << "I have found " << *nprimes << " primes in " << static_cast<double>(elapsed.count())/1000 << " seconds in the set " << nstart << ".\n";
    delete [] P;
}
