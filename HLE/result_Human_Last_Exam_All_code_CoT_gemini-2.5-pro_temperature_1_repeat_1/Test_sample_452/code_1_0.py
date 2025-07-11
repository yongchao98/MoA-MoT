import math

# The problem asks for the constant b in the asymptotic formula
# C(n) ~ b * n^(3/2) for the expected cover-and-return time on a
# random n-vertex tree.
# Based on established results in probability theory, primarily from
# David Aldous, the constant b is given by sqrt(pi / 2).

# Define the constants
pi_val = math.pi
constant_b = math.sqrt(pi_val / 2)

# The problem asks to output the final equation with each number.
print("The constant b is the leading coefficient in the asymptotic formula for the expected cover time of a random tree.")
print("The exact value of this constant is given by the expression:")
print("b = sqrt(pi / 2)")
print("")
print("Let's calculate the numerical value:")
print(f"Using pi = {pi_val}")
print(f"b = sqrt({pi_val} / 2)")
print(f"b = {constant_b}")
