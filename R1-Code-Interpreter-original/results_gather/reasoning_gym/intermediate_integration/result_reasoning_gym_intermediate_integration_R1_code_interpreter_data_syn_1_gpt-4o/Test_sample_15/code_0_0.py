from sympy import symbols, atan, integrate, ln

# Define the variable and function
x = symbols('x')
f = -3 * atan(2 * x)

# Compute the antiderivative
antiderivative = integrate(f, x)

# Print the result
print(antiderivative)