import sympy

# Let the connective constant be mu.
# Based on known results for 4-regular lattices in statistical mechanics,
# the most plausible value for the connective constant is 1 + sqrt(2).
# We will find the minimal polynomial for this value.
x = sympy.Symbol('x')
# The value mu = 1 + sqrt(2) is a root of the polynomial.
# Let x = 1 + sqrt(2)
# x - 1 = sqrt(2)
# (x - 1)^2 = 2
# x^2 - 2x + 1 = 2
# x^2 - 2x - 1 = 0
# The polynomial is P(x) = x^2 - 2x - 1.
# This polynomial is monic, has rational coefficients, and is irreducible over Q
# because its other root is 1 - sqrt(2), which is irrational.
# Therefore, it is the minimal polynomial.

# The coefficients of the polynomial x^2 - 2x - 1 = 0 are a=1, b=-2, c=-1.
a = 1
b = -2
c = -1

# As requested, printing out each number in the final equation.
print(f"The minimal polynomial is P(x) = a*x^2 + b*x + c = 0")
print(f"The final equation is: {a}*x^2 + ({b})*x + ({c}) = 0")
