from sympy import symbols, asin, sqrt, integrate

# Define the variable
x = symbols('x')

# Define the function to integrate
f = -2 * asin(x)

# Compute the antiderivative
antiderivative = integrate(f, x)

# Print the result
print(antiderivative)