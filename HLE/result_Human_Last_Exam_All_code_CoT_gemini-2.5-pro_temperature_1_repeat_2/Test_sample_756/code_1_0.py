from fractions import Fraction

# This script calculates and prints the maximum value of |b| + |c|
# based on the derived optimal polynomial f(x) = ax^2 + bx + c.

# The coefficients of the optimal polynomial are:
a = Fraction(-8, 9)
b = Fraction(8, 9)
c = Fraction(7, 9)

# The expression to maximize is |b| + |c|.
abs_b = abs(b)
abs_c = abs(c)
max_value = abs_b + abs_c

print(f"An optimal polynomial that maximizes |b| + |c| is f(x) = ({a})x^2 + ({b})x + ({c}).")
print("For this polynomial, we calculate the value of |b| + |c|.")
print("\nThe final equation is:")
print(f"|{b}| + |{c}| = {abs_b} + {abs_c} = {max_value}")