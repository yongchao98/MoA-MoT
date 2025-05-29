from sympy import symbols, integrate

# Define the variable
x = symbols('x')

# Define the function
function = -5/2*x**4 + 8*x**3 + 18/7*x**2

# Calculate the antiderivative
antiderivative = integrate(function, x)

# Print the result
print(antiderivative)