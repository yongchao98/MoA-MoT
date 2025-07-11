import numpy as np

def find_smallest_prime_factor(k):
    """Finds the smallest prime factor of an integer k."""
    k = abs(int(k))
    if k < 2:
        return None # No prime factors
    if k % 2 == 0:
        return 2
    for i in range(3, int(k**0.5) + 1, 2):
        if k % i == 0:
            return i
    return k

# The figure eight knot has 4 arcs and 4 crossings.
# Let the arcs be a, b, c, d. A standard diagram gives the following
# equations (mod n) at the crossings:
# 1: a + d = 2b  =>  a - 2b + d = 0
# 2: b + c = 2d  =>  b + c - 2d = 0
# 3: c + d = 2a  => -2a + c + d = 0
# 4: a + b = 2c  =>  a + b - 2c = 0
#
# This gives the coloring matrix M:
#      a   b   c   d
M = np.array([
    [ 1, -2,  0,  1],  # eq 1
    [ 0,  1,  1, -2],  # eq 2
    [-2,  0,  1,  1],  # eq 3
    [ 1,  1, -2,  0]   # eq 4
])

# The knot determinant is the absolute value of the determinant of any
# (n-1)x(n-1) minor of the coloring matrix. We'll remove the last row and column.
M_minor = M[:-1, :-1]

# Calculate the knot determinant.
knot_determinant = np.linalg.det(M_minor)

# The smallest n for a non-trivial coloring is the smallest prime factor
# of the knot determinant.
n = find_smallest_prime_factor(knot_determinant)

print(f"The coloring matrix of the figure eight knot is:\n{M}")
print(f"\nA minor of this matrix is:\n{M_minor}")
print(f"\nThe determinant of the minor is: {knot_determinant:.0f}")
print(f"The knot determinant is the absolute value: {abs(int(knot_determinant))}")
print(f"\nThe smallest n for a non-trivial coloring is the smallest prime factor of the determinant.")
print(f"Therefore, n = {n}.")
print("\nThis means the smallest algebraic structure, Z_n, has 5 elements: {0, 1, 2, 3, 4}.")

# We can find an example of a valid 5-coloring. Let a=0 and b=1.
# From the equations, we can deduce a coloring (a, b, c, d) = (0, 1, 3, 2).
# Let's verify this with one of the crossing equations: a + b = 2c (mod 5).
a, b, c = 0, 1, 3
lhs = a + b
rhs = 2 * c
print(f"\nTo satisfy the final equation requirement, we show an example.")
print(f"An example of a crossing equation is: a + b = 2 * c (mod {n})")
print(f"Using a valid coloring with a={a}, b={b}, c={c}:")
print(f"The final equation is: {a} + {b} = 2 * {c} (mod {n})")
print(f"Which simplifies to: {lhs} = {rhs} (mod {n})")
print(f"This holds, since {rhs} % {n} is {rhs % n}, which equals {lhs}.")

print(f"\nThe number of elements in the smallest algebraic structure that allows coloring the figure eight knot is {n}.")
