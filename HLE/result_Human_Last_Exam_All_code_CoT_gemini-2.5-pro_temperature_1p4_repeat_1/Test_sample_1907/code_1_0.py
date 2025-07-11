import math

def get_nth_prime(n):
    """
    Calculates the n-th prime number.
    It estimates an upper bound for the n-th prime, uses a sieve to find all primes
    up to that bound, and then returns the n-th one.
    """
    if n <= 0:
        return None
    if n == 1:
        return 2

    # Estimate an upper bound for the n-th prime.
    # For n >= 6, p_n < n * (ln(n) + ln(ln(n))).
    # We use a slightly more generous bound to be safe.
    if n < 6:
        # A simple bound for small n
        limit = 15
    else:
        log_n = math.log(n)
        # Add a factor to ln(ln(n)) to ensure the bound is sufficient
        limit = math.ceil(n * (log_n + 1.15 * math.log(log_n)))

    # Sieve of Eratosthenes to find primes up to the limit
    primes_bool = [True] * (limit + 1)
    primes_bool[0] = primes_bool[1] = False
    for i in range(2, int(math.sqrt(limit)) + 1):
        if primes_bool[i]:
            for multiple in range(i * i, limit + 1, i):
                primes_bool[multiple] = False

    # Create a list of prime numbers from the boolean array
    prime_numbers = [i for i, is_p in enumerate(primes_bool) if is_p]

    if len(prime_numbers) < n:
        # Fallback if the estimate was too low, though unlikely with the chosen formula
        print(f"Warning: Sieve limit {limit} was too small for n={n}. Recalculating with a larger limit.")
        return get_nth_prime_robust(n)

    return prime_numbers[n - 1]

def get_nth_prime_robust(n):
    """A slower but guaranteed method if the sieve limit fails."""
    primes = [2]
    num = 3
    while len(primes) < n:
        is_prime = True
        for p in primes:
            if p * p > num:
                break
            if num % p == 0:
                is_prime = False
                break
        if is_prime:
            primes.append(num)
        num += 2
    return primes[-1]


# Let x_k be the first term of P^(k)
# The sequence is defined by x_1 = 2 and x_k = p_{x_{k-1}} for k >= 2

# Start with x_1
current_term = 2
print(f"The 1st term in P^(1) is {current_term}")

# Iterate from k=2 to 11 to find x_11
for k in range(2, 12):
    index = current_term
    current_term = get_nth_prime(index)
    print(f"The 1st term in P^({k}) is p_{index} = {current_term}")

print(f"\nThe final answer for the 1st term in P^(11) is {current_term}")
<<<9878201>>>