from sympy import symbols, integrate

# Define the variable and the function
X = symbols('X')
function = 2*X**9 - 30*X**4 + 3

# Calculate the antiderivative
antiderivative = integrate(function, X)

# Print the result
print(antiderivative)