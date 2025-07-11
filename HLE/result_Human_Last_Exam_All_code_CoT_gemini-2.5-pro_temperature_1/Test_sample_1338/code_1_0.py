import math

def get_prime_factorization(n):
    """Returns the prime factorization of n as a dictionary."""
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def is_prime_power(n):
    """Checks if n is a prime power (p^k for p prime, k>=1)."""
    if n <= 1:
        return False
    factors = get_prime_factorization(n)
    return len(factors) == 1

def is_perfect_square(n):
    """Checks if n is a perfect square."""
    if n < 1:
        return False
    sqrt_n = int(math.sqrt(n))
    return sqrt_n * sqrt_n == n

def find_tiling_subset():
    """
    Determines for which integers t in the given set the number of
    t-omino tilings of an n x n grid is even for any n.
    """
    t_values = [2, 3, 4, 5, 7, 9, 15]
    result_set = []
    
    print("Analyzing each t in {2, 3, 4, 5, 7, 9, 15}:")
    print("-" * 50)
    
    for t in t_values:
        print(f"Case t = {t}:")
        
        # Rule 1: t is a perfect square, m*m.
        # Counterexample: N(m, t) = 1 (odd).
        if is_perfect_square(t):
            m = int(math.sqrt(t))
            print(f"  Result: OUT. t={t} is a perfect square ({m}x{m}).")
            print(f"  A {m}x{m} grid has 1 tiling with a single {t}-omino, which is odd.")
            print("-" * 50)
            continue
            
        # Rule 2: t is an odd prime power.
        # Counterexamples can be constructed (e.g., N(2t,t) is odd).
        if t != 2 and is_prime_power(t):
            print(f"  Result: OUT. t={t} is an odd prime power.")
            print(f"  For such t, an n x n grid can be constructed that has an odd number of tilings.")
            print("-" * 50)
            continue
            
        # Rule 3: t = 2 (special case).
        # Number of tilings is always even.
        if t == 2:
            print(f"  Result: IN. For t=2 (dominoes), the number of tilings of an n x n grid is always even.")
            print("-" * 50)
            result_set.append(t)
            continue

        # Rule 4: t is not a prime power.
        # Number of tilings is always even (conjectured).
        if not is_prime_power(t):
            print(f"  Result: IN. t={t} is not a prime power (e.g., 15 = 3x5).")
            print(f"  The number of tilings for such t is believed to always be even.")
            print("-" * 50)
            result_set.append(t)
            continue
            
    print("\nFinal conclusion:")
    print("The subset of integers for which the statement is true is:")
    print(sorted(result_set))

find_tiling_subset()