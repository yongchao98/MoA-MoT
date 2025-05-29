from sympy import symbols, integrate, asin, sqrt

# Define the variable and function
x = symbols('x')
f = -3 * asin(2 * x)

# Compute the antiderivative
antiderivative = integrate(f, x)

# Output the result
print(antiderivative)