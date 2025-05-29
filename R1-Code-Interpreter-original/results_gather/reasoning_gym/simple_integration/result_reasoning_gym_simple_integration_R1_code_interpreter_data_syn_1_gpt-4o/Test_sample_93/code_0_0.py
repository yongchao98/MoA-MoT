from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = (25/3)*x**4 - (21/2)*x**2

# Calculate the antiderivative
antiderivative = integrate(function, x)

# Print the result
print(antiderivative)