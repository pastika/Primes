#include <ctime>
#include <cmath>
#include <iostream>
#include <cstdint>

//Type of variable to bit pack
typedef std::uint32_t T;
//Number of numbers to check (must be multiple of sizeof(T)*16)
const std::uint64_t NEVTS = 1000000000;

template< std::uint64_t N> class LOG2
{
public:
    static const std::uint64_t RESULT = 1 + LOG2<N/2>::RESULT;
};

template<> class LOG2< 2 >
{
public:
    static const std::uint64_t RESULT = 1;
};

//These values are automatically calculated
const std::uint64_t SA = (std::uint64_t)sqrt((double)NEVTS)/2;


template<std::uint64_t NEVTS, typename T>
class BitArray
{
private:
    //These values are automatically calculated
    static constexpr std::uint64_t N = sizeof(T)*8;
    static constexpr std::uint64_t A = NEVTS/(N*2);
    static constexpr std::uint64_t L2N = LOG2<N>::RESULT;

    T P[A];

public:
    BitArray(bool d)
    {
        if(d) for(std::uint64_t i = 0; i < A; i++) P[i] = ~(T)(0);
        else  for(std::uint64_t i = 0; i < A; i++) P[i] = (T)(0);
    }
    
    inline bool operator[](std::uint64_t i) const
    {
        return bool(P[i >> L2N] & (T(1) << (i & (N - 1))));
    }

    template<bool v = false>
    inline void assign(const std::uint64_t i)
    {
        if(v) P[i >> L2N] |=   T(1) << (i & (N - 1));
        else  P[i >> L2N] &= ~(T(1) << (i & (N - 1)));
    }

    inline std::uint64_t count()
    {
        std::uint64_t count = 0 ;
        for(auto n : P)
        {
            while (n)
            {
                ++count;
                n &= (n - 1);     
            }
        }
        return count;
    }
};

int main()
{
    clock_t t0, t1;
    unsigned int numPrimes = 0, k, i;
        
    t0 = clock();
    static BitArray<NEVTS, T> P(true);
    for(k = 1; k <= SA; ++k) if(P[k]) for(i = 2*k*(k + 1); i < NEVTS/2; i += (2*k + 1)) P.assign<false>(i);

    numPrimes = P.count();
    t1 = clock();

    std::cout << "I have found " << numPrimes + 1 << " primes in " << (double)(t1 - t0)/CLOCKS_PER_SEC << " seconds.\n";
}

