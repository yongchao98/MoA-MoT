import sys

def find_prime_solution():
    """
    Searches for a prime number n such that for some k > 2,
    n is the m-th prime number where m = (n+1)/k, and m is also prime.
    """
    # Set a search limit for the primes.
    # A limit of 1,000,000 allows testing for n up to 1,000,000 and m up to pi(1,000,000) ~ 78,498.
    # This is a reasonably large search space.
    LIMIT = 1000000

    # Step 1: Generate primes up to LIMIT using a Sieve of Eratosthenes.
    is_prime = [True] * (LIMIT + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(LIMIT**0.5) + 1):
        if is_prime[i]:
            for multiple in range(i*i, LIMIT + 1, i):
                is_prime[multiple] = False
    
    primes = [i for i, is_p in enumerate(is_prime) if is_p]
    
    num_primes = len(primes)
    
    # Step 2: Iterate through primes m and check the conditions.
    # The iterator `p` will be our prime candidate for `m`.
    for m in primes:
        # We need to find the m-th prime, which requires our list to have at least m elements.
        # This check ensures we don't go out of bounds.
        if m > num_primes:
            break
            
        # n is the m-th prime number. Our list is 0-indexed, so we access index m-1.
        n = primes[m - 1]
        
        # Check if n+1 is divisible by m
        if (n + 1) % m == 0:
            k = (n + 1) // m
            
            # Check if k > 2
            if k > 2:
                # Solution found. Print the details and exit.
                print(f"Solution found for k > 2:")
                print(f"The prime n is {n}.")
                print(f"The index m = (n+1)/k is {m}, which is prime.")
                print(f"The value for k is {k}.")
                print(f"The final equation is {n} + 1 = {k} * {m}.")
                print("\nAnswer:")
                print(f"{k}, {n}")
                return

    # If the loop completes, no solution was found within the search limit.
    print("Does not exist")

find_prime_solution()