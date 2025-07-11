def find_special_prime():
    """
    Finds the largest prime p < 1,000,000 of the form p = 4u+1,
    where u is a prime of the form u = 4v+1, and v is also a prime.
    """
    limit = 1000000
    
    # Step 1: Use Sieve of Eratosthenes to find all primes up to the limit.
    is_prime = [True] * limit
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            for multiple in range(i * i, limit, i):
                is_prime[multiple] = False
    
    # Step 2: Iterate downwards from the limit to find the largest p.
    # The condition p=4u+1 and u=4v+1 implies p = 4(4v+1)+1 = 16v+5.
    # So we only need to check primes p where p % 16 == 5.
    for p in range(limit - 1, 1, -1):
        # Step 3: Check if p satisfies the conditions.
        # Condition a: p is prime.
        # Condition b: p = 16v + 5, which means (p - 5) must be divisible by 16.
        if is_prime[p] and (p - 5) % 16 == 0:
            v = (p - 5) // 16
            u = 4 * v + 1 # or u = (p-1)//4
            
            # Since p is within the limit, u and v will also be.
            # No need for boundary checks on u and v.
            
            # Step 3c: Check if u and v are also prime.
            if is_prime[u] and is_prime[v]:
                # Found the largest p that satisfies all conditions.
                # Print the numbers p, u, and v.
                print(str(p) + ":" + str(u) + ":" + str(v))
                return

# Execute the function to find and print the result.
find_special_prime()