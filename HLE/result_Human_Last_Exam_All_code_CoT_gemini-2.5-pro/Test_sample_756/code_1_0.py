import numpy as np

# Based on the derivation, the maximum value of |b| + |c| is 1.5.
# This maximum is achieved for polynomials like P(x) = -0.5*x^2 - x + 0.5.
# Let's define the coefficients of this polynomial.
a = -0.5
b = -1.0
c = 0.5

# We can verify that for this polynomial, |P(x)| <= 1 for all x in [-1, 1].
# P'(x) = -x - 1. For x in [-1, 1], P'(x) <= 0, so P(x) is decreasing.
# P(-1) = -0.5*(-1)^2 - (-1) + 0.5 = -0.5 + 1 + 0.5 = 1.
# P(1) = -0.5*(1)^2 - 1 + 0.5 = -0.5 - 1 + 0.5 = -1.
# Since the function is monotonic on [-1, 1], its values are between -1 and 1.
# The condition is satisfied.

# Now, we calculate the value of |b| + |c|.
abs_b = abs(b)
abs_c = abs(c)
max_value = abs_b + abs_c

print("The problem is to find the maximum value of |b| + |c| for a quadratic P(x) = ax^2+bx+c")
print("such that |P(x)| <= 1 for all x in [-1, 1].")
print("\nAn optimal polynomial that maximizes this expression is P(x) = -0.5x^2 - x + 0.5.")
print(f"For this polynomial, the coefficients are a = {a}, b = {b}, c = {c}.")

print("\nThe value of the expression |b| + |c| is calculated as follows:")
print(f"|b| + |c| = |{b}| + |{c}| = {abs_b} + {abs_c} = {max_value}")

print(f"\nThe maximum value of |b| + |c| is {max_value}.")
