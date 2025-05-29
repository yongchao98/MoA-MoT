from sympy import symbols, integrate

# Define the variable
X = symbols('X')

# Define the function
function = -18*X**8 - 6*X

# Calculate the antiderivative
antiderivative = integrate(function, X)

# Print the result
print(antiderivative)