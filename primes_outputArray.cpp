#include <cmath>
#include <ctime>
#include <iostream>
#include <vector>
#include <cstdint>

constexpr unsigned int NEVTS = 1000000000;

std::vector<int> getPrimes(const unsigned int NEVTS)
{
    const int A = NEVTS / 2;
    const int SA = (int)sqrt((double)NEVTS) / 2;

    bool *sieve = new bool[A];
    std::vector<int> primes;
    //preload 2 because we removed all even numbers
    primes.push_back(2);

    for (unsigned int i = 0; i < A; ++i) sieve[i] = true;
    for (unsigned int primeToCheck = 1; primeToCheck <= SA; ++primeToCheck)
    {
        if (sieve[primeToCheck])
        {
            for (unsigned int numToMask = 2 * primeToCheck*(primeToCheck + 1); numToMask < A; numToMask += (2 * primeToCheck + 1))
            {
                sieve[numToMask] = false;
            }
        }
    }


    for (unsigned int i = 1; i < A; i++)
    {
        if (sieve[i])
        {
            primes.push_back(i * 2 + 1);
        }
    }

    delete[] sieve;

    return primes;
}

int main()
{
    clock_t t0, t1;
    t0 = clock();
    std::vector<int> primes = getPrimes(NEVTS);
    t1 = clock();

    std::cout << "I have found " << primes.size() << " primes in " << (double)(t1 - t0) / CLOCKS_PER_SEC << " seconds.\n";
}
