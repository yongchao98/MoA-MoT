from sympy import symbols, integrate, sin, cos

# Define the variable and the function
x = symbols('x')
f = -6 * sin(2*x + 2) * cos(2*x + 2)

# Calculate the antiderivative
antiderivative = integrate(f, x)

# Print the result
print(antiderivative)