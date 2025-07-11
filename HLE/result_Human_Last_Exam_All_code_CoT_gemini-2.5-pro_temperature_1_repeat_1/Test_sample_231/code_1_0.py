import math

def solve_prime_problem():
    """
    Finds the smallest integer k > 2 and prime n such that n is the m-th prime,
    where m = (n+1)/k is also a prime number.
    This is equivalent to finding the smallest k > 2 and prime m such that p_m = k*m - 1.
    """

    # We need a list of primes to find p_m. The list needs to be large enough.
    # p_m is roughly m*log(m). For the search, m can be up to e^k.
    # For k=10, m can be around 22000, and p_m around 250000.
    # A sieve limit of 300,000 should be sufficient to find a solution for small k.
    LIMIT = 300000
    
    # Step 1: Generate primes using Sieve of Eratosthenes
    primes_bool = [True] * (LIMIT + 1)
    primes_bool[0] = primes_bool[1] = False
    for i in range(2, int(math.sqrt(LIMIT)) + 1):
        if primes_bool[i]:
            for multiple in range(i * i, LIMIT + 1, i):
                primes_bool[multiple] = False
    
    primes_list = [i for i, is_prime in enumerate(primes_bool) if is_prime]
    primes_set = set(primes_list)

    # Step 2: Search for k and m
    # Outer loop for k starting from 3
    for k in range(3, 100): # A reasonable upper bound for k
        # Inner loop for prime m
        for m in primes_list:
            # The m-th prime is primes_list[m-1].
            # This requires our list to have at least m elements.
            if m > len(primes_list):
                # This would mean our LIMIT is too small.
                # For the expected solution, it's not an issue.
                break

            # The value of the m-th prime, p_m
            p_m = primes_list[m - 1]
            
            # The potential value of n based on the formula
            n = k * m - 1

            # Optimization: if p_m gets larger than n, it will likely stay larger
            # for this k, so we can move to the next k.
            if p_m > n:
                break
            
            # Check if we found a solution
            if p_m == n:
                # All conditions are met by construction:
                # k > 2 is from the loop.
                # m is a prime number from primes_list.
                # n = p_m, so n is a prime number.
                # (n+1)/k = (p_m+1)/k = (km-1+1)/k = m, which is prime.
                
                print("Found a solution!")
                print(f"The equation is n = p_((n+1)/k)")
                print(f"The smallest values are k = {k} and n = {n}.")
                m_val = (n + 1) // k
                print(f"Verification:")
                print(f"k = {k}")
                print(f"n = {n}")
                print(f"m = (n+1)/k = ({n}+1)/{k} = {m_val}")
                print(f"m = {m_val} is a prime number.")
                print(f"The {m_val}-th prime number is indeed {p_m}.")
                print("\nFinal Answer:")
                print(f"{k}, {n}")
                return

    print("Does not exist within the searched range.")

# Run the solver
solve_prime_problem()