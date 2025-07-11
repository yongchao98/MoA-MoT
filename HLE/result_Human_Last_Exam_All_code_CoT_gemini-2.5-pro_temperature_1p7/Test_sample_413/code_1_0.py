from fractions import Fraction

# Based on our analysis, the polynomial is f(x) = x^3 - (3/8)x^2 - (5/8)x.
# Let's define the coefficients a, b, c from f(x) = x^3 + ax^2 + bx + c.
a = Fraction(-3, 8)
b = Fraction(-5, 8)
c = Fraction(0)

# We want to compute f(3).
x = 3

# Calculate each term of the polynomial for x=3.
term1 = x**3
term2 = a * (x**2)
term3 = b * x
term4 = c

# The sum of the terms gives the final result.
result = term1 + term2 + term3 + term4

# Output the computation, showing each number in the final equation.
print(f"The determined polynomial is f(x) = x^3 + ({a})x^2 + ({b})x + {c}")
print(f"To compute f(3), we substitute x=3 into the expression:")
print(f"f(3) = 3^3 + ({a}) * 3^2 + ({b}) * 3 + {c}")
print(f"The values of the terms are:")
print(f"3^3 = {term1}")
print(f"({a}) * 3^2 = {term2}")
print(f"({b}) * 3 = {term3}")
print(f"c = {term4}")
print(f"The final sum is:")
print(f"f(3) = {term1} + ({term2}) + ({term3}) = {result}")
