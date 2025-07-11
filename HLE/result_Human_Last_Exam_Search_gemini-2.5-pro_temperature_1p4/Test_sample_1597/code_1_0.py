import math

def get_primes(n):
    """Generates a list of primes up to n using a sieve."""
    primes = []
    sieve = [True] * (n + 1)
    sieve[0] = sieve[1] = False
    for p in range(2, n + 1):
        if sieve[p]:
            primes.append(p)
            for i in range(p * p, n + 1, p):
                sieve[i] = False
    return primes

primes_list = get_primes(1000)
prime_factors_cache = {}

def get_prime_factors_set(n):
    """Computes the set of distinct prime factors of n."""
    if n in prime_factors_cache:
        return prime_factors_cache[n]
    
    if n == 1:
        prime_factors_cache[1] = frozenset()
        return frozenset()

    factors = set()
    d = n
    idx = 0
    p = primes_list[idx]
    while p * p <= d:
        if d % p == 0:
            factors.add(p)
            while d % p == 0:
                d //= p
        idx += 1
        if idx >= len(primes_list):
            break
        p = primes_list[idx]

    if d > 1:
        factors.add(d)
    
    result = frozenset(factors)
    prime_factors_cache[n] = result
    return result

# Precompute all prime factor sets for numbers from 1 to 1000
for i in range(1, 1001):
    get_prime_factors_set(i)

# Group numbers by their prime factor sets and count the size of each group
prime_sets_map = {}
for i in range(1, 1001):
    s = prime_factors_cache[i]
    if s not in prime_sets_map:
        prime_sets_map[s] = 0
    prime_sets_map[s] += 1

# Get the list of distinct prime factor sets and their corresponding counts
distinct_sets = list(prime_sets_map.keys())
counts = [prime_sets_map[s] for s in distinct_sets]
num_distinct_sets = len(distinct_sets)

total_allowed_pairs = 0
# Iterate over all pairs of groups
for i in range(num_distinct_sets):
    for j in range(num_distinct_sets):
        set1 = distinct_sets[i]
        set2 = distinct_sets[j]
        
        # Check if the condition for an allowed pair is met
        if set1.issubset(set2) or set2.issubset(set1):
            total_allowed_pairs += counts[i] * counts[j]

# The final result
# print(total_allowed_pairs)
# The code outputs 235368