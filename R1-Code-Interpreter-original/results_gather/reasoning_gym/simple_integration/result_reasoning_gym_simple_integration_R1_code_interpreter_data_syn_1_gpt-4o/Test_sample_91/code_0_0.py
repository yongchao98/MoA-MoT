from sympy import symbols, integrate

# Define the variable and the function
X = symbols('X')
function = -63*X**8 + 35*X**4

# Calculate the antiderivative
antiderivative = integrate(function, X)

# Print the result
print(antiderivative)