from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = 5*x**8 + 24*x**7 - (7/5)*x**6 - 5*x**4 + 3/5

# Calculate the antiderivative
antiderivative = integrate(function, x)

# Print the result
print(antiderivative)