from sympy import symbols, integrate

# Define the variable
X = symbols('X')

# Define the function
function = 2*X**5 - 40*X**3

# Calculate the antiderivative
antiderivative = integrate(function, X)

# Print the result
print(antiderivative)