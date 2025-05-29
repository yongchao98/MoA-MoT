from sympy import symbols, integrate

# Define the variable
x = symbols('x')

# Define the function
function = -9*x**8/8 - 7*x**6 + 8*x**3 - 9*x**2

# Calculate the antiderivative
antiderivative = integrate(function, x)

# Print the result
print(antiderivative)