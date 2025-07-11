import math

def get_prime_factorization(n):
    """Computes the prime factorization of n."""
    factors = {}
    temp = n
    d = 2
    while d * d <= temp:
        if temp % d == 0:
            count = 0
            while temp % d == 0:
                count += 1
                temp //= d
            factors[d] = count
        d += 1
    if temp > 1:
        factors[temp] = 1
    return factors

def gcd_list(numbers):
    """Computes the greatest common divisor of a list of numbers."""
    if not numbers:
        return 1
    result = numbers[0]
    for i in range(1, len(numbers)):
        result = math.gcd(result, numbers[i])
    return result

def get_root(n):
    """
    Computes the integer root of n.
    The root of n is r such that n = r^k and r is not a perfect power.
    This is found by dividing all exponents in the prime factorization of n by their GCD.
    """
    if n == 1:
        return 1
    
    factors = get_prime_factorization(n)
    exponents = list(factors.values())
    
    g = gcd_list(exponents)
    
    if g == 1:
        return n
        
    root = 1
    for p, e in factors.items():
        root *= p**(e // g)
    return root

def solve():
    """
    Calculates the total number of allowed pairs (a,b) with 1 <= a,b <= 1000.
    """
    LIMIT = 1000

    # Case 1: Pairs where a=1 or b=1.
    # There are LIMIT pairs (1,b), LIMIT pairs (a,1). The pair (1,1) is in both sets.
    count_pairs_with_1 = LIMIT + LIMIT - 1
    
    # Case 2: Pairs where a > 1 and b > 1.
    # These pairs are allowed if root(a) == root(b).
    # We group numbers by their root and sum the squares of the group sizes.
    root_counts = {}
    for n in range(2, LIMIT + 1):
        r = get_root(n)
        root_counts[r] = root_counts.get(r, 0) + 1
        
    squared_counts = []
    # Sort by root for a deterministic output order
    for r in sorted(root_counts.keys()):
        c = root_counts[r]
        squared_counts.append(c * c)
    
    total_count = count_pairs_with_1 + sum(squared_counts)
    
    # Print the equation with all its components, as requested.
    print(f"{count_pairs_with_1} + ", end="")
    print(" + ".join(map(str, squared_counts)), end="")
    print(f" = {total_count}")

solve()