def find_prime_solution():
    """
    Searches for a prime number n and an integer k>2 satisfying the given conditions.

    Conditions:
    1. n is a prime number.
    2. k is an integer > 2.
    3. m = (n+1)/k is a prime number.
    4. n is the m-th prime number (P_m).
    """
    # Step 1: Set a reasonable search limit for the Sieve.
    # The number of primes up to x is pi(x) ~ x/ln(x).
    # To test prime indices m, we need P_m. We need pi(SIEVE_LIMIT) > m.
    # pi(400000) is ~33860. So we can test prime indices m up to 33860.
    SIEVE_LIMIT = 400000

    # Step 2: Generate primes using Sieve of Eratosthenes.
    is_prime_sieve = [True] * (SIEVE_LIMIT + 1)
    is_prime_sieve[0] = is_prime_sieve[1] = False
    for p in range(2, int(SIEVE_LIMIT**0.5) + 1):
        if is_prime_sieve[p]:
            for multiple in range(p * p, SIEVE_LIMIT + 1, p):
                is_prime_sieve[multiple] = False

    primes = [p for p, is_p in enumerate(is_prime_sieve) if is_p]
    num_generated_primes = len(primes)
    
    solutions = []

    # Step 3: Iterate through possible indices 'm'.
    # The m-th prime has index m, which must be prime itself.
    # The index `m` (1-based) cannot be larger than the number of primes we've generated.
    # We loop m_idx (0-based) from 0 upwards. The m-th prime is primes[m-1].
    for m_idx in range(num_generated_primes):
        m = m_idx + 1  # m is the 1-based index we are testing.

        # Condition: The index 'm' must be prime.
        if m < 2 or not is_prime_sieve[m]:
            continue

        # At this point, m is a prime number.
        # Condition: n is the m-th prime number.
        n = primes[m_idx]

        # From m = (n+1)/k, we get k = (n+1)/m.
        # Check if k is an integer and k > 2.
        if (n + 1) % m == 0:
            k = (n + 1) // m
            if k > 2:
                # We found a valid solution (k, n).
                solutions.append((k, n))

    # Step 4: Report the best solution or "Does not exist".
    if not solutions:
        print("Does not exist")
    else:
        # Sort by k first, then by n to find the smallest pair.
        solutions.sort()
        best_k, best_n = solutions[0]
        m_val = (best_n + 1) // best_k
        print(f"Smallest solution found: k={best_k}, n={best_n}")
        print(f"Verification:")
        print(f"n = {best_n} is prime.")
        print(f"The index m = (n+1)/k = ({best_n}+1)/{best_k} = {m_val}, which is prime.")
        print(f"The {m_val}-th prime number is indeed {best_n}.")


if __name__ == '__main__':
    find_prime_solution()