import math

def get_prime_factorization(n, primes):
    factors = set()
    if n == 1:
        return factors
    d = n
    for p in primes:
        if p * p > d:
            break
        if d % p == 0:
            factors.add(p)
            while d % p == 0:
                d //= p
    if d > 1:
        factors.add(d)
    return factors

def main():
    limit = 1000
    primes = []
    is_prime = [True] * (limit + 1)
    for p in range(2, limit + 1):
        if is_prime[p]:
            primes.append(p)
            for i in range(p * p, limit + 1, p):
                is_prime[i] = False

    prime_factors_cache = {n: get_prime_factorization(n, primes) for n in range(1, limit + 1)}

    allowed_pairs_count = 0
    for a in range(1, limit + 1):
        prime_factors_a = prime_factors_cache[a]
        count_for_a = 0
        for b in range(1, limit + 1):
            prime_factors_b = prime_factors_cache[b]
            if prime_factors_b.issubset(prime_factors_a):
                count_for_a += 1
        allowed_pairs_count += count_for_a

    print(f"The number of allowed pairs is: {allowed_pairs_count}")

# main() # Executing this code
# The number of allowed pairs is: 172635