#include <RooRealVar.h>
#include "RooFit/MultiProcess/Config.h"

#include <sstream>

#include "gtest/gtest.h"

bool isPrime(long number)
{
	for (long a = 2; a <= sqrt(number); a++)
	{
		if (number % a == 0)
		{
			return false;
		}
	}
	return true;
}

int compute_primes(int amount, long min_value)
{
	int total_primes(0);
	for (long currentNum = min_value;; currentNum++)
	{
		if (isPrime(currentNum))
		{
			//std::cout << endl << currentNum << " ";
			total_primes++;
			if (total_primes == amount) return total_primes;
		}
	}
	//return total_primes;
	return 0;
}

TEST(Test, test)
{
    RooFit::MultiProcess::Config::callgrind_zero();
    std::cout << "printing message " << std::endl;
    pid_t pid = getpid();

    compute_primes(1e5, 1e6);

    printf("pid: %lu\n", pid);
    RooFit::MultiProcess::Config::callgrind_dump();

   EXPECT_TRUE(0 == 0);
}