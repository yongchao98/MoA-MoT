from sympy import symbols, integrate, sqrt

# Define the variable and the function
x = symbols('x')
f = 20 * sqrt(-10*x - 6)

# Calculate the antiderivative
antiderivative = integrate(f, x)

# Print the result
print(antiderivative)