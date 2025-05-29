from sympy import symbols, integrate

# Define the variable and the function
X = symbols('X')
function = 25*X**(9/2) - 21*X**6 + 24*X**2 - 35/4

# Calculate the antiderivative
antiderivative = integrate(function, X)

# Print the result
print(antiderivative)