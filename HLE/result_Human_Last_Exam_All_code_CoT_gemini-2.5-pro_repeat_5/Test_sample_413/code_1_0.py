from fractions import Fraction

# From the derivation, the coefficients of the polynomial f(x) = x^3 + ax^2 + bx + c are:
a = Fraction(-3, 8)
b = Fraction(-5, 8)
c = Fraction(0)

# The value of x to be evaluated
x = 3

# Calculate f(3)
f3 = x**3 + a * x**2 + b * x + c

# Print the final equation with each number as requested
print(f"The function is determined to be f(x) = x^3 + ({a})x^2 + ({b})x.")
print(f"To compute f(3), we substitute x=3 into the function:")
print(f"f(3) = 3^3 + ({a}) * 3^2 + ({b}) * 3")

# Showing the calculation steps
term1 = x**3
term2 = a * x**2
term3 = b * x

print(f"f(3) = {term1} + ({term2}) + ({term3})")
print(f"f(3) = {term1} - {abs(term2)} - {abs(term3)}")
print(f"f(3) = 27 - 27/8 - 15/8")
print(f"f(3) = 27 - (27 + 15)/8")
print(f"f(3) = 27 - 42/8")
print(f"f(3) = 27 - 21/4")
print(f"f(3) = 108/4 - 21/4")
print(f"f(3) = {f3.numerator}/{f3.denominator}")
print(f"The final exact value of f(3) is: {f3}")