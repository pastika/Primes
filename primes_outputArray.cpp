#include <cmath>
#include <ctime>
#include <iostream>
#include <vector>

const int NEVTS = 1000000000;
const int A = NEVTS/2;
const int SA = (int)sqrt((double)NEVTS)/2;

std::vector<int> getPrimes(int NEVTS, int A, int SA)
{
	unsigned int k, i;
	bool *P = new bool[A];
	std::vector<int> primes;
	primes.push_back(2);

	for (i = 0; i < A; ++i) P[i] = true;
	for (k = 1; k <= SA; ++k)
	{
		if (P[k])
		{
			//primes.push_back(k * 2 + 1);
			for (i = 2 * k*(k + 1); i < A; i += (2 * k + 1))
			{
				P[i] = false;
			}
		}
	}
	

	for (i = 1; i < A; i++)
	{
		if (P[i])
		{
			primes.push_back(i * 2 + 1);
		}
	}

	delete[] P;

	return primes;
}

int main()
{
	clock_t t0, t1;
	t0 = clock();
	std::vector<int> primes = getPrimes(NEVTS, A, SA);
	t1 = clock();
	
	std::cout << "I have found " << primes.size() << " primes in " << (double)(t1 - t0)/CLOCKS_PER_SEC << " seconds.\n";

	//for (int i = 0; i < primes.size(); ++i)
	//{
	//	std::cout << primes[i] << std::endl;
	//}
}