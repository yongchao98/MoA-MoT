import math

def get_distinct_prime_factors(n):
    """
    Returns a set of distinct prime factors of a number n.
    """
    factors = set()
    # Check for number of 2s
    while n % 2 == 0:
        factors.add(2)
        n = n // 2
    # n must be odd at this point. Check for factors from 3 up to sqrt(n)
    for i in range(3, int(math.sqrt(n)) + 1, 2):
        while n % i == 0:
            factors.add(i)
            n = n // i
    # If n is a prime number greater than 2
    if n > 2:
        factors.add(n)
    return factors

def is_prime_power(t):
    """
    Checks if a number t is a prime power (p^k for p prime, k>=1).
    """
    if t <= 1:
        return False
    
    distinct_factors = get_distinct_prime_factors(t)
    
    return len(distinct_factors) == 1

def solve():
    """
    Finds the subset of integers for which the number of t-omino tilings
    of an n x n grid is always even. This is conjectured to be the set of
    prime powers.
    """
    numbers = [2, 3, 4, 5, 7, 9, 15]
    result_subset = []
    for t in numbers:
        if is_prime_power(t):
            result_subset.append(t)
    
    print("The subset of integers is:", result_subset)

solve()
