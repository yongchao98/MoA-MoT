from sympy import symbols, integrate

# Define the variable
X = symbols('X')

# Define the function
f = 16*X**7 - 5*X**4/3 + 10*X

# Calculate the antiderivative
antiderivative = integrate(f, X)

# Print the result
print(antiderivative)