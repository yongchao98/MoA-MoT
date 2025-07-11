def solve_prime_puzzle():
    """
    Searches for a prime number n and an integer k > 2 such that
    m = (n+1)/k is prime and n is the m-th prime number.
    """
    # Set a limit for the Sieve of Eratosthenes. A larger limit allows for a more
    # extensive search. A limit of 200,000 is a reasonable starting point.
    # p_m is approximately m*ln(m). To test m up to ~17000, we need primes up to
    # 17000 * ln(17000) ~= 17000 * 9.74 ~= 165580. So 200,000 is sufficient.
    LIMIT = 200000

    # 1. Generate primes using Sieve of Eratosthenes
    is_prime = [True] * (LIMIT + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(LIMIT**0.5) + 1):
        if is_prime[i]:
            for multiple in range(i * i, LIMIT + 1, i):
                is_prime[multiple] = False
    
    primes = [i for i, is_p in enumerate(is_prime) if is_p]
    num_primes = len(primes)

    # 2. Iterate through candidate primes 'm'
    # The value 'm' must be prime, so we iterate through our generated list.
    for m_candidate in primes:
        m = m_candidate
        
        # We need to find the m-th prime, which is primes[m-1].
        # We must ensure our list is long enough for this index.
        if m - 1 >= num_primes:
            # Our sieve limit is too small to find the m-th prime.
            # We have exhausted the searchable range.
            break

        # 3. Define n as the m-th prime
        n = primes[m - 1]
        
        # 4. Check if (n+1) is divisible by m
        if (n + 1) % m == 0:
            # 5. Calculate k
            k = (n + 1) // m
            
            # 6. Check if k > 2
            if k > 2:
                # A solution is found. Since we iterate m in increasing order,
                # this will be the solution with the smallest n and likely smallest k.
                print(f"{k}, {n}")
                return

    # 7. If the loop completes, no solution was found in the searched range.
    print("Does not exist")

solve_prime_puzzle()