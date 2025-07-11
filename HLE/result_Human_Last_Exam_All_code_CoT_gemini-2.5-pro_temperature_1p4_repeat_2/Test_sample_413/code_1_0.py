from fractions import Fraction

# Coefficients determined from the analysis
a = Fraction(-3, 2)
b = Fraction(-19, 16)
c = Fraction(27, 32)

# Define the function f(x)
def f(x_val):
    x = Fraction(x_val)
    return x**3 + a * x**2 + b * x + c

# The value to compute is f(3)
x_in = 3
result = f(x_in)

# Print the final equation with fractions
x = Fraction(x_in)
term1 = x**3
term2 = a * x**2
term3 = b * x
term4 = c

print(f"Based on the analysis, the function is f(x) = x^3 + ({a})x^2 + ({b})x + ({c}).")
print(f"We need to compute f(3).")
print(f"f(3) = (3)^3 + ({a})*(3)^2 + ({b})*(3) + ({c})")
print(f"f(3) = {term1} + ({term2}) + ({term3}) + {term4}")
print(f"f(3) = {term1} {term2} {term3} + {term4}")
print(f"The exact value of f(3) is: {result}")
