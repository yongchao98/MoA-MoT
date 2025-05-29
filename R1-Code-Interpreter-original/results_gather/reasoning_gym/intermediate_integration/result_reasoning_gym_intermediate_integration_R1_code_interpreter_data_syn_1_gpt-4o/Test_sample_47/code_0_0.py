from sympy import symbols, integrate, sqrt

# Define the variable and the function
x = symbols('x')
f = 4 * sqrt(4*x + 8)

# Calculate the antiderivative
antiderivative = integrate(f, x)

# Print the result
print(antiderivative)