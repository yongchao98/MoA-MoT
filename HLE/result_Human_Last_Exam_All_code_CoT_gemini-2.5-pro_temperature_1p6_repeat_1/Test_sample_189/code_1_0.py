def find_prime_chain():
    """
    This function finds the largest prime p < 1,000,000 such that:
    1. p = 4u + 1, where u is prime.
    2. u = 4v + 1, where v is prime.
    """
    
    limit = 999999
    
    # Step 1: Use a Sieve of Eratosthenes to find all primes up to the limit.
    # is_prime[i] will be True if i is a prime number, False otherwise.
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    for p in range(2, int(limit**0.5) + 1):
        if is_prime[p]:
            for multiple in range(p * p, limit + 1, p):
                is_prime[multiple] = False

    largest_p_found = -1
    result_u = -1
    result_v = -1

    # Step 2: Iterate through possible primes 'v'.
    # From p = 16v + 5, if p < limit, then v < limit / 16.
    for v in range(2, limit // 16 + 2):
        # We only proceed if v is a prime number.
        if is_prime[v]:
            # Step 3: Calculate u and check if it's prime.
            u = 4 * v + 1
            
            # If u is out of bounds, we can stop.
            if u > limit:
                break
            
            if is_prime[u]:
                # Step 4: Calculate p and check if it's prime.
                p = 4 * u + 1
                
                # If p is out of bounds, continue to the next v.
                if p > limit:
                    continue
                
                if is_prime[p]:
                    # A valid chain is found. We update the largest p.
                    if p > largest_p_found:
                        largest_p_found = p
                        result_u = u
                        result_v = v

    # Step 5: Print the final result in the format p:u:v.
    if largest_p_found != -1:
        print(f"{largest_p_found}:{result_u}:{result_v}")
    else:
        print("No such prime chain was found within the specified limit.")

# Execute the function to find and print the result.
find_prime_chain()
<<<958553:239638:59909>>>