import math

# The problem asks for the closed expression of the infinite product
# P = product_{n=0 to inf} (1 - exp(-(2n+1)pi))
# Based on the theory of elliptic functions and Jacobi theta functions,
# this product evaluates to 2^(-1/8).

# Here we will demonstrate this value.
result_val = 2**(-1/8)

# The closed expression can be written in several ways.
# We will print the final equation.

numerator = 1
denominator_base = 2
denominator_root = 8

print(f"The closed-form expression for the infinite product is 2^(-1/8).")
print(f"This can be written as:")
print(f"P = {denominator_base}^(-1/{denominator_root})")
print(f"Which is equal to:")
print(f"P = {numerator} / {denominator_base}^(1/{denominator_root})")
print(f"In radical form:")
print(f"P = {numerator} / \u221B{denominator_base}")
print(f"The numerical value is approximately: {result_val}")
