#include <cmath>
#include <vector>
#include <algorithm>
#include <ctime>
#include <iostream>
#include <cstdint>

const std::uint64_t NEVTS = 1000000000ul;
const std::uint64_t A = NEVTS/2;
const std::uint64_t SA = (std::uint64_t)sqrt((double)NEVTS)/2;

int main()
{
    clock_t t0, t1;
    std::uint64_t numPrimes = 0, k, i;
        
    t0 = clock();
    std::vector<bool> P(A, true);
    for(k = 1; k <= SA; ++k) if(P[k]) for(i = 2*k*(k+1); i < A; i += (2*k + 1)) P[i] = false;

    numPrimes = std::count(P.begin(), P.end(), true);
    t1 = clock();
        
    std::cout << "I have found " << numPrimes + 1 << " primes in " << (double)(t1 - t0)/CLOCKS_PER_SEC << " seconds.\n";
}
