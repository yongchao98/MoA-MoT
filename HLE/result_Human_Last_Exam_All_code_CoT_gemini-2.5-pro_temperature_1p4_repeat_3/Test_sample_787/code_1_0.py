import math

def get_primes_up_to(n):
    """Generates a list of primes up to n using a sieve."""
    primes = []
    is_prime = [True] * (n + 1)
    is_prime[0] = is_prime[1] = False
    for p in range(2, n + 1):
        if is_prime[p]:
            primes.append(p)
            for i in range(p * p, n + 1, p):
                is_prime[i] = False
    return primes

def valuation(n, q):
    """Calculates the q-adic valuation of n."""
    if n == 0:
        return float('inf')
    if q <= 1:
        return 0
    count = 0
    while n != 0 and n % q == 0:
        count += 1
        n //= q
    return count

def P(x):
    """Evaluates the polynomial P(X)."""
    # Use integer arithmetic to maintain precision
    res = 1
    for k in range(5):
        term = pow(x, 5) - pow(x, k)
        res *= term
    return res

# The analysis shows the only prime factors of the limit L can be 2, 3, and 5.
# We find the exponent of each prime factor in L by finding the minimum
# valuation of P(p) for that prime factor, across various primes p.
prime_factors_of_L = [2, 3, 5]
exponents = {}

# We test with primes up to 100, which is sufficient to find the minimums.
test_primes = get_primes_up_to(100)

for q in prime_factors_of_L:
    min_vq = float('inf')
    # We search for min v_q(P(p)) for p != q.
    for p in test_primes:
        if p == q:
            continue
        
        val = valuation(P(p), q)
        if val < min_vq:
            min_vq = val
    exponents[q] = min_vq

v_2_L = exponents[2]
v_3_L = exponents[3]
v_5_L = exponents[5]

limit_g_n = pow(2, v_2_L) * pow(3, v_3_L) * pow(5, v_5_L)

print("The limit L is the product of prime powers q^(v_q(L)).")
print(f"The exponent for q=2 is v_2(L) = {v_2_L}")
print(f"The exponent for q=3 is v_3(L) = {v_3_L}")
print(f"The exponent for q=5 is v_5(L) = {v_5_L}")
print("\nThe limit is calculated as:")
print(f"L = 2^{v_2_L} * 3^{v_3_L} * 5^{v_5_L} = {pow(2,v_2_L)} * {pow(3,v_3_L)} * {pow(5,v_5_L)} = {limit_g_n}")
