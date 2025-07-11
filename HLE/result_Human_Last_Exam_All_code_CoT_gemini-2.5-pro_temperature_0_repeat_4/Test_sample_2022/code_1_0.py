import math

def combinations(n, k):
    """Calculates the binomial coefficient C(n, k)"""
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

def solve():
    """
    Calculates the value of F(n) for the given primes.
    The value n modulo p is 6. So we calculate F(6).
    S(6) = sum_{k=0 to 3} C(6-2k, k) * (-1/8)^k
    The terms are:
    k=0: C(6,0) * (-1/8)^0 = 1
    k=1: C(4,1) * (-1/8)^1 = -4/8
    k=2: C(2,2) * (-1/8)^2 = 1/64
    k=3: C(0,3) * (-1/8)^3 = 0
    So, F(6) = (1 - 4 * 8^-1 + 1 * 64^-1) mod p
    """
    primes = [80039, 80077]
    results = []

    for p in primes:
        # We need to calculate S(6) mod p
        n_val = 6
        
        # Term for k=0
        # C(6,0) * (-1/8)^0 = 1
        term0_val = combinations(n_val - 2*0, 0)
        term0 = term0_val
        
        # Term for k=1
        # C(4,1) * (-1/8)^1 = -4 * (8^-1)
        term1_val = combinations(n_val - 2*1, 1)
        inv_8 = pow(8, p - 2, p)
        term1 = (p - (term1_val * inv_8) % p) % p

        # Term for k=2
        # C(2,2) * (-1/8)^2 = 1 * (64^-1)
        term2_val = combinations(n_val - 2*2, 2)
        inv_64 = pow(64, p - 2, p)
        term2 = (term2_val * inv_64) % p
        
        # Term for k=3 is 0 since C(0,3) = 0
        
        final_result = (term0 + term1 + term2) % p
        results.append(final_result)

    print(f"The final equation is S(6) mod p = (C(6,0) - C(4,1) * 8^-1 + C(2,2) * 64^-1) mod p")
    print(f"For p = {primes[0]}: ({combinations(6,0)} - {combinations(4,1)} * {pow(8, primes[0]-2, primes[0])} + {combinations(2,2)} * {pow(64, primes[0]-2, primes[0])}) mod {primes[0]} = {results[0]}")
    print(f"For p = {primes[1]}: ({combinations(6,0)} - {combinations(4,1)} * {pow(8, primes[1]-2, primes[1])} + {combinations(2,2)} * {pow(64, primes[1]-2, primes[1])}) mod {primes[1]} = {results[1]}")
    print(f"Final answers separated by a comma: {results[0]},{results[1]}")

solve()