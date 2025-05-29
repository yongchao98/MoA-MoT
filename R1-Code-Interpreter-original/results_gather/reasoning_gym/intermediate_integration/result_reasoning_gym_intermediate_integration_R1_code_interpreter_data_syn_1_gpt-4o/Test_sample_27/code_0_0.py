from sympy import symbols, integrate, sqrt

# Define the variable and the function
x = symbols('x')
f = 9 * sqrt(3*x - 1)

# Calculate the antiderivative
antiderivative = integrate(f, x)

# Print the result
print(antiderivative)