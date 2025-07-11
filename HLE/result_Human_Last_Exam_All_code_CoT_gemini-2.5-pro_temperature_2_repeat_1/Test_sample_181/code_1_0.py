import math

# This script calculates foo(7) by determining the mathematical equivalent of the function.
# The function foo(n) counts the number of distinct multiplication tables on n elements
# that satisfy the axioms of a group. This is equivalent to finding the number of
# labeled groups of order n.

# For n = 7, since 7 is a prime number, any group of order 7 is isomorphic
# to the cyclic group Z_7.

# The number of such labeled groups is given by the formula: n! / |Aut(G)|,
# where G is the group structure (Z_7 in this case) and |Aut(G)| is the
# order of its automorphism group.

# The order of the automorphism group of Z_n is phi(n), Euler's totient function.
# For a prime number p, phi(p) = p - 1.

n = 7

# Step 1: Calculate n!
factorial_n = math.factorial(n)

# Step 2: Calculate |Aut(Z_7)| = phi(7)
num_automorphisms = n - 1

# Step 3: Calculate the final result
result = factorial_n // num_automorphisms

# Print the final equation with each number
print(f"foo(7) = {n}! / phi({n})")
print(f"foo(7) = {factorial_n} / {num_automorphisms}")
print(f"foo(7) = {result}")
