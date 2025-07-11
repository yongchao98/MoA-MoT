from fractions import Fraction
from math import factorial

# Use a cache for memoization to speed up Bernoulli number calculation
_bernoulli_cache = {}

def combinations(n, k):
    """Calculates the binomial coefficient C(n, k)"""
    if k < 0 or k > n:
        return 0
    return factorial(n) // (factorial(k) * factorial(n - k))

def bernoulli(n):
    """
    Calculates the n-th Bernoulli number using the recurrence relation.
    B_n = - sum_{k=0}^{n-1} (C(n+1, k) * B_k) / (n+1)
    """
    if n in _bernoulli_cache:
        return _bernoulli_cache[n]
    
    if n == 0:
        return Fraction(1)
    
    sum_val = Fraction(0)
    for k in range(n):
        sum_val += combinations(n + 1, k) * bernoulli(k)
    
    result = -sum_val / (n + 1)
    _bernoulli_cache[n] = result
    return result

def zeta_rational_part(n):
    """
    Calculates the rational part of zeta(n) for an even integer n.
    zeta(n) / pi^n = (-1)^(k+1) * B_n * 2^(n-1) / n!  (where n=2k)
    """
    if n % 2 != 0 or n == 0:
        raise ValueError("Input must be a positive even integer.")
        
    k = n // 2
    b_n = bernoulli(n)
    
    sign = 1 if (k + 1) % 2 == 0 else -1
    
    numerator = sign * b_n.numerator * (2**(n - 1))
    denominator = b_n.denominator * factorial(n)
    
    return Fraction(numerator, denominator)

def solve_sum():
    """
    Calculates the sum by computing zeta(6)*zeta(8)/zeta(24).
    """
    # The sum is equal to zeta(6) * zeta(8) / zeta(24)
    # The power of pi is 6 + 8 - 24 = -10.
    
    # We calculate the rational part of each zeta function value.
    z6_rat = zeta_rational_part(6)
    z8_rat = zeta_rational_part(8)
    z24_rat = zeta_rational_part(24)
    
    # The final rational coefficient is z6_rat * z8_rat / z24_rat
    final_rational_coeff = z6_rat * z8_rat / z24_rat
    
    num = final_rational_coeff.numerator
    den = final_rational_coeff.denominator
    pi_power = -10
    
    # The problem asks to output each number in the final equation.
    print(f"{num} / {den} * pi^({pi_power})")

solve_sum()