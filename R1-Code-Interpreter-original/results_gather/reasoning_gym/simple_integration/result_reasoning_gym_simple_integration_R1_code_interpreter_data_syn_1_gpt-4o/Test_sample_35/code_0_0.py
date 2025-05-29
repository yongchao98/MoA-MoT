from sympy import symbols, integrate

# Define the variable
x = symbols('x')

# Define the function to integrate
function = 2/3 - 35/9 * x**4

# Calculate the antiderivative
antiderivative = integrate(function, x)

# Print the result
print(antiderivative)