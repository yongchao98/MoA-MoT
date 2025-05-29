from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = -54*x**5 + 5*x**4 + 18*x**2 + 9

# Calculate the antiderivative
antiderivative = integrate(function, x)

# Print the result
print(antiderivative)