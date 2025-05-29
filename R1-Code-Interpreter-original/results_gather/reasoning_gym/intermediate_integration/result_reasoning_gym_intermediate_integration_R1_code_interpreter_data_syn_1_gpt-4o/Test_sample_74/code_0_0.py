from sympy import symbols, exp, integrate

# Define the variable and the function
x = symbols('x')
function = exp(-x - 2)

# Calculate the antiderivative
antiderivative = integrate(function, x)

# Print the result
print(antiderivative)