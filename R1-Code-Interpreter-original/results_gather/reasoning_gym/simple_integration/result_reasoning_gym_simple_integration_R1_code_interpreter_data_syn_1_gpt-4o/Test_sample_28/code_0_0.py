from sympy import symbols, integrate

# Define the variable
X = symbols('X')

# Define the function
function = -20*X**3 + 21*X**2 + 2*X/7 + 7

# Compute the antiderivative
antiderivative = integrate(function, X)

# Print the result
print(antiderivative)