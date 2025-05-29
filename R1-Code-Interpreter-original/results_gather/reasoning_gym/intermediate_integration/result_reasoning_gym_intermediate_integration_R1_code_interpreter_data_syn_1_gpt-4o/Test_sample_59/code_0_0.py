from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = -2 * (3*x - 7)**4

# Calculate the antiderivative
antiderivative = integrate(function, x)

# Print the result
print(antiderivative)