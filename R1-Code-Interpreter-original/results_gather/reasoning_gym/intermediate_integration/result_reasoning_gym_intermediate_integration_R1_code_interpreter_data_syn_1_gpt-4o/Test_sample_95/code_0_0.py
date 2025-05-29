from sympy import symbols, integrate

# Define the variable
x = symbols('x')

# Define the function to integrate
function = 2 * (x - 1)**2

# Calculate the antiderivative
antiderivative = integrate(function, x)

# Print the result
print(antiderivative)