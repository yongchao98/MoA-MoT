import math
from collections import defaultdict

def combinations(n, k):
    """Computes the binomial coefficient C(n, k)"""
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    if k > n // 2:
        k = n - k
    
    res = 1
    for i in range(k):
        res = res * (n - i) // (i + 1)
    return res

def get_prime_factorization(num):
    """Returns a dictionary of prime factors of a number"""
    factors = defaultdict(int)
    d = 2
    temp = abs(num)
    while d * d <= temp:
        while temp % d == 0:
            factors[d] += 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] += 1
    return factors

def solve():
    """
    Calculates the maximum order among all elements of the Smith group of S.
    """
    n = 50
    k = 5

    # Step 1: Compute eigenvalues of the adjacency matrix A
    P_j2 = []
    for j in range(k + 1):
        val = (combinations(5 - 0, 2 - 0) * combinations(5 - j, 0) * combinations(n - k - j + 0, 0) -
               combinations(5 - 1, 2 - 1) * combinations(5 - j, 1) * combinations(n - k - j + 1, 1) +
               combinations(5 - 2, 2 - 2) * combinations(5 - j, 2) * combinations(n - k - j + 2, 2))
        P_j2.append(val)

    # Step 2: Compute eigenvalues of the Seidel matrix S
    v = combinations(n, k)
    s = []
    # Eigenvalue for the all-ones vector
    s_0 = v - 1 - 2 * P_j2[0]
    s.append(s_0)
    # Other eigenvalues
    for j in range(1, k + 1):
        s_j = -1 - 2 * P_j2[j]
        s.append(s_j)

    # Step 3 & 4: Find prime factors and maximum exponents for each prime
    prime_max_exponents = defaultdict(int)
    for eig in s:
        if eig == 0: continue
        factors = get_prime_factorization(eig)
        for prime, exponent in factors.items():
            if exponent > prime_max_exponents[prime]:
                prime_max_exponents[prime] = exponent
    
    # Step 5: Calculate the final maximum order
    max_order = 1
    factors_str = []
    
    # Sort primes for deterministic output
    sorted_primes = sorted(prime_max_exponents.keys())

    for prime in sorted_primes:
        exponent = prime_max_exponents[prime]
        term = prime ** exponent
        max_order *= term
        factors_str.append(f"{prime}^{exponent}")
    
    print("The maximum order is the product of the following prime powers:")
    print(" * ".join(factors_str))
    print(f"\nThis evaluates to:\n{max_order}")

solve()