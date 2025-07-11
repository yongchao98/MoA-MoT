import sympy

# The particular solution derived is x(t) = -8 / (4*t^2 - 8*t + 1)
# We need to find the value of x(1).

t = 1

# Numerator of the expression
numerator = -8

# Calculate the terms in the denominator
term1 = 4 * (t**2)
term2 = -8 * t
term3 = 1

# Calculate the denominator
denominator = term1 + term2 + term3

# Calculate the final result
# Using sympy.Rational to keep the result as a fraction
result = sympy.Rational(numerator, denominator)

# Print the calculation steps
print(f"The particular solution is x(t) = -8 / (4*t^2 - 8*t + 1).")
print(f"We substitute t = {t} to find x({t}):")
print(f"x({t}) = {numerator} / (4*({t})^2 - 8*({t}) + 1)")
print(f"x({t}) = {numerator} / ({term1} - {abs(term2)} + {term3})")
print(f"x({t}) = {numerator} / ({denominator})")
print(f"x({t}) = {result}")