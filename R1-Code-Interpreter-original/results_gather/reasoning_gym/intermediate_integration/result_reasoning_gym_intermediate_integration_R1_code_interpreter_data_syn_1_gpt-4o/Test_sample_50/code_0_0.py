from sympy import symbols, integrate, sqrt

# Define the variable and the function
x = symbols('x')
f = -3*sqrt(5 - 3*x)

# Calculate the antiderivative
antiderivative = integrate(f, x)

# Print the result
print(antiderivative)