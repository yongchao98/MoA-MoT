from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = -25*x**4/7 - 4*x**3 + 10*x/3

# Calculate the antiderivative
antiderivative = integrate(function, x)

# Print the result
print(antiderivative)