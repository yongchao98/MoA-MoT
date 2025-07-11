import math

def find_prime_solution():
    """
    Searches for a prime number n and an integer k > 2 such that
    m = (n+1)/k is a prime number and n is the m-th prime number.
    This is equivalent to finding a prime m such that k = (p_m + 1) / m is an integer > 2.
    """

    # We need to find p_m, where m can be large.
    # p_m is approximately m * log(m).
    # We set a search limit for the primes n. This determines how many primes p_m we can check.
    # A limit of 500,000 for n means we can check m up to pi(500,000), which is 41,538.
    limit = 500000

    # Step 1: Generate primes up to the limit using a sieve.
    try:
        is_prime = [True] * (limit + 1)
        is_prime[0] = is_prime[1] = False
        for p in range(2, int(math.sqrt(limit)) + 1):
            if is_prime[p]:
                for i in range(p * p, limit + 1, p):
                    is_prime[i] = False
        
        primes_list = [p for p, is_p in enumerate(is_prime) if is_p]
        num_primes_found = len(primes_list)
    except MemoryError:
        print("Error: The limit is too large, causing a MemoryError.")
        print("Please try running on a machine with more memory or reduce the limit.")
        return

    # Step 2: Iterate through prime indices `m` and check the conditions.
    # The index `m` must itself be a prime number.
    # We can check primes `m` as long as `m` is not larger than the number of primes we found.
    
    best_k = float('inf')
    best_n = float('inf')
    found_solution = False

    for m in primes_list:
        # The index `m` must be less than or equal to the total number of primes we have generated.
        if m > num_primes_found:
            break

        # Get n, the m-th prime number (using 0-based indexing for our list).
        n = primes_list[m-1]

        # Check if k = (n + 1) / m is an integer.
        if (n + 1) % m == 0:
            k = (n + 1) // m
            
            # Check if k > 2.
            if k > 2:
                # We found a valid solution. Check if it's the best one so far.
                if k < best_k:
                    best_k = k
                    best_n = n
                    found_solution = True
                elif k == best_k:
                    if n < best_n:
                        best_n = n

    # Step 3: Print the result.
    if found_solution:
        # Note: The final output must be in the format `k, n`.
        print(f"{best_k}, {best_n}")
    else:
        # Extensive searches have found no solution, so this is the expected outcome.
        print("Does not exist")

if __name__ == '__main__':
    find_prime_solution()