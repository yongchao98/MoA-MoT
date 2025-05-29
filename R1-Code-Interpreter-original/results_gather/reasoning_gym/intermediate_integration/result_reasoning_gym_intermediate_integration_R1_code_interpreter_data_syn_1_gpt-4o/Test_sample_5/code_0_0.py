from sympy import symbols, exp, integrate

# Define the variable and the function
x = symbols('x')
f = 24 * exp(8 - 8*x)

# Calculate the antiderivative
antiderivative = integrate(f, x)

# Print the result
print(antiderivative)