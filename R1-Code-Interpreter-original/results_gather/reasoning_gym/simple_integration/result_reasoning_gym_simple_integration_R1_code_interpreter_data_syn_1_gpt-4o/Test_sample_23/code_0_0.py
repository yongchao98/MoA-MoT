from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = -16*x**7 + 10*x

# Calculate the antiderivative
antiderivative = integrate(function, x)

# Print the result
print(antiderivative)