from sympy import symbols, integrate

# Define the variable and the function
X = symbols('X')
function = 64*X**7 + 30*X**5 + 5

# Calculate the antiderivative
antiderivative = integrate(function, X)

# Print the result
print(antiderivative)