from sympy import symbols, integrate

# Define the variable and the function
X = symbols('X')
function = -5 * X**9

# Calculate the antiderivative
antiderivative = integrate(function, X)

# Print the result
print(antiderivative)