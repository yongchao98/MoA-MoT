from sympy import symbols, integrate

# Define the variable
X = symbols('X')

# Define the function
f = 80*X**9 - 90*X**8 + 2*X/3 + 1

# Compute the antiderivative
antiderivative = integrate(f, X)

# Print the result
print(antiderivative)