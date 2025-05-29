from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = x**8 - 24*x**7 - 35*x**6 - 16*x

# Calculate the antiderivative
antiderivative = integrate(function, x)

# Print the result
print(antiderivative)