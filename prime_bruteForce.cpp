#include <iostream>
#include <cmath>
#include <ctime>
#include <vector>

int main()
{
	using namespace std;

	clock_t t0, t1;
	t0 = clock();

	vector<int> primes;
	primes.push_back(2);
	vector<int>::const_iterator j;

	for(int i = 3; i < 1000000000; i+=2)
	{
		//const int sqrti = (int)sqrt(double(i));
		for(j = primes.begin(); (i%(*j)) && j != primes.end(); ++j)	
		{
			if((*j)*(*j) > i)
			{
				primes.push_back(i);
				break;
			}
		}
		/*for(const auto& j : primes)	
		{
			if(!(i%j)) break;
			else if(j > sqrti)
			{
				primes.push_back(i);
				break;
			}
		}*/
	}
	t1 = clock();

	std::cout << "I have found " << primes.size() + 1 << " primes in " << (double)(t1 - t0)/CLOCKS_PER_SEC << " seconds.\n";
}