import math

# A global cache to store primes so we don't have to re-sieve repeatedly.
_primes_cache = []

def get_nth_prime(n):
    """
    Returns the n-th prime number.
    Uses a cached list of primes and extends it with a sieve if necessary.
    """
    global _primes_cache
    if not isinstance(n, int) or n <= 0:
        raise ValueError("n must be a positive integer.")

    if n <= len(_primes_cache):
        return _primes_cache[n - 1]

    # Estimate the upper limit for the sieve to find the n-th prime.
    # The Prime Number Theorem provides approximations. A common upper bound is:
    # p_n < n * (ln(n) + ln(ln(n))) for n >= 6.
    if n < 6:
        # Small n, a small limit is sufficient.
        limit = 15
    else:
        log_n = math.log(n)
        log_log_n = math.log(log_n)
        limit = math.ceil(n * (log_n + log_log_n))
    
    # Run the sieve until we have enough primes.
    # In case the estimate is too low, we increase it and retry.
    while len(_primes_cache) < n:
        # The sieve creates a boolean array up to the limit.
        is_prime = [True] * (limit + 1)
        is_prime[0] = is_prime[1] = False
        
        # Sieve of Eratosthenes
        for i in range(2, int(math.sqrt(limit)) + 1):
            if is_prime[i]:
                for multiple in range(i * i, limit + 1, i):
                    is_prime[multiple] = False
        
        # Rebuild the prime cache from the sieve result.
        new_primes = [i for i, p_bool in enumerate(is_prime) if p_bool]
        
        # If the sieve was large enough, we are done. Otherwise, increase the limit.
        if len(new_primes) >= n:
            _primes_cache = new_primes
        else:
            # Increase limit by 20% and try again
            limit = int(limit * 1.2)


    return _primes_cache[n - 1]

def solve():
    """
    Solves the problem by iteratively calculating the first term of P^(k).
    """
    print("Let a_k be the 1st term in P^(k).")
    print("The recurrence relation is a_k = p(a_{k-1}), where p(n) is the n-th prime.")
    print("-" * 40)
    
    # a_1 = 2, the first prime number.
    sequence_terms = [2]
    current_term = 2
    
    print(f"a_1 = {current_term}")

    # Calculate a_2 through a_11
    for k in range(2, 12):
        prev_term = current_term
        try:
            current_term = get_nth_prime(prev_term)
            sequence_terms.append(current_term)
            print(f"a_{k} = p(a_{k-1}) = p({prev_term}) = {current_term}")
        except Exception as e:
            print(f"An error occurred while calculating a_{k}: {e}")
            return

    final_answer = sequence_terms[-1]
    print("-" * 40)
    print(f"The 1st term in P^(11) is a_11.")
    print(f"Final Answer: {final_answer}")
    
solve()
<<<8969711>>>