import math

# The optimal position for the vertex of the parabola is found to be x_v = 1/2.
x_v = 1/2

# Calculate the coefficients a, b, and c for the polynomial P(x) = ax^2 + bx + c
# that maximizes |b|+|c| under the given constraints.
# a = -2 / (1 + x_v)**2
b = (4 * x_v) / (1 + x_v)**2
c = (1 + 2 * x_v - x_v**2) / (1 + x_v)**2

# The maximum value is the sum of the absolute values of b and c.
max_value = abs(b) + abs(c)

# We print the calculation step by step.
# Using fractions for precision in the explanation
from fractions import Fraction
b_frac = Fraction(b).limit_denominator()
c_frac = Fraction(c).limit_denominator()
max_val_frac = Fraction(max_value).limit_denominator()

print(f"The optimal polynomial has coefficients b = {b_frac} and c = {c_frac}.")
print("The maximum value of |b| + |c| is calculated as:")
print(f"|{b_frac}| + |{c_frac}| = {b_frac} + {c_frac} = {max_val_frac}")
print(f"In floating point, the value is {max_value}")
