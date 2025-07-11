import math

# Based on the analysis, the problem asks for a well-known constant in optimization theory.
# The smallest value z is the worst-case ratio of the SDP relaxation value to the
# true value of a quadratic program, where the supremum is taken over all
# positive semidefinite "cost" matrices.
# This constant is proven to be pi / 2.

# We will now print the equation and its numerical value.
equation_string = "z = pi / 2"
pi_value = math.pi
z_value = pi_value / 2

print("The smallest value of z is determined by the equation:")
print(equation_string)
print(f"z = {pi_value} / 2")
print(f"z = {z_value}")