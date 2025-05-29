from sympy import symbols, integrate

# Define the variable
X = symbols('X')

# Define the function
function = (-20/3)*X**9 - 36*X**7 - (7/8)*X**6 + 12*X**3 - 12*X

# Calculate the antiderivative
antiderivative = integrate(function, X)

# Print the result
print(antiderivative)