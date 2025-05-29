from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = -35*x**6/4 + 48*x**5 - 45*x**4 + 6

# Calculate the antiderivative
antiderivative = integrate(function, x)

# Print the result
print(antiderivative)