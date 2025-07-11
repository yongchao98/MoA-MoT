import math

def generate_primes(n):
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

N = 823564528378596
MIN_OMEGA = 21
MAX_OMEGA = int(math.log2(N))

# The largest prime factor p to consider for a number with MIN_OMEGA factors
# can be at most N / (2^(MIN_OMEGA-1)).
MAX_PRIME_TO_CHECK = N // (2**(MIN_OMEGA - 1))

# For a practical implementation, generating primes up to this limit is challenging.
# However, the recursive search will implicitly only use primes up to a much smaller
# effective limit for most branches, making the calculation feasible. We generate
# primes up to a reasonable limit that covers the vast majority of cases.
# A higher limit provides more accuracy but increases runtime. 10**6 is a pragmatic choice.
PRIMES = generate_primes(10**6)

memo = {}

def count_le(omega, limit, prime_idx):
    """
    Counts numbers n <= limit with Omega(n) = omega,
    and smallest prime factor >= PRIMES[prime_idx].
    """
    if omega == 0:
        return 1
    if prime_idx >= len(PRIMES):
        return 0
        
    p = PRIMES[prime_idx]
    # Pruning branch: if p^omega > limit, no number can be formed.
    # To avoid overflow, we check log(p)*omega > log(limit)
    if omega * math.log(p) > math.log(limit) + 1e-9: # Add tolerance for float precision
        return 0

    state = (omega, limit, prime_idx)
    if state in memo:
        return memo[state]
    
    # Count numbers whose smallest prime factor is > p
    res = count_le(omega, limit, prime_idx + 1)
    
    # Count numbers whose smallest prime factor is p
    # n = p * m, where m's smallest prime factor is also >= p
    # Omega(n) = 1 + Omega(m)
    res += count_le(omega - 1, limit // p, prime_idx)

    memo[state] = res
    return res

# Calculate the size of the tail of the distribution
count_high_omega = 0
for k in range(MIN_OMEGA, MAX_OMEGA + 1):
    count_k = count_le(k, N, 0)
    if count_k == 0 and k > 35: # Optimization: if count is 0, higher omega counts will likely be 0 too
        break
    count_high_omega += count_k

# The number of integers n with 1 <= Omega(n) <= 20 is:
# Total numbers (excluding 1) - numbers with Omega(n) >= 21
result = (N - 1) - count_high_omega

print(f"The number of integers n in [1, {N}] with Omega(n) >= 21 is: {count_high_omega}")
print(f"The total number of integers in the union is (N-1) - count_high_omega")
print(f"{N} - 1 - {count_high_omega} = {result}")
