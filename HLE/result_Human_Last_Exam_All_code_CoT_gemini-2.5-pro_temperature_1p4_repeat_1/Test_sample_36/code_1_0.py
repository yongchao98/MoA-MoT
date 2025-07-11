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

# The number of unique sets is determined by the properties of the group with the
# largest exponent. For Abelian groups of order 18, the possible exponents are 18 (for C_18)
# and 6 (for C_2 x C_3 x C_3). The largest exponent is 18.

# The problem reduces to counting the non-empty order ideals in the divisor lattice of 18.
# The prime factorization of 18 is 2^1 * 3^2.
# The exponents of the prime factors are 1 and 2.

# The divisor lattice of a number n = p1^a1 * p2^a2 * ... is isomorphic to the product 
# of chains C_{a1+1} x C_{a2+1} x ...
# For 18, the chains are of sizes m = 1+1=2 and n = 2+1=3.

# The number of order ideals in a product of two chains of size m and n is C(m+n, m).
m = 1 + 1
n = 2 + 1
num_ideals = combinations(m + n, m)
num_unique_sets = num_ideals - 1

print("The number of unique sets is found by the following calculation:")
print("1. Find the maximum exponent of Abelian groups of order 18. This is 18 (from C_18).")
print("2. Count the non-empty order ideals of the divisor lattice of 18.")
print("3. The prime factorization of 18 is 2^1 * 3^2.")
print("4. The divisor lattice is isomorphic to a product of chains of size (1+1) and (2+1).")
print(f"5. The number of ideals is given by the binomial coefficient C((1+1)+(2+1), 1+1) = C({m+n}, {m}).")
print("\nFinal Equation:")
print(f"C({m}, {n}) = C({m+n}, {m}) = {num_ideals}")
print(f"The number of unique non-empty sets is one less than the total number of ideals.")
print(f"Result = {num_ideals} - 1 = {num_unique_sets}")
<<<9>>>