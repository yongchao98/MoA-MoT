from sympy import symbols, integrate, sqrt

# Define the variable and the function
x = symbols('x')
f = -15 * sqrt(7 - 5 * x)

# Calculate the antiderivative
antiderivative = integrate(f, x)

# Print the result
print(antiderivative)