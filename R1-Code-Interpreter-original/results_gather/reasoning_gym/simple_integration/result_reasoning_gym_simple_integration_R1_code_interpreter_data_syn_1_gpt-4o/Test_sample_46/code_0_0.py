from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = (56/5)*x**6 - 18*x**5 + 7/6

# Calculate the antiderivative
antiderivative = integrate(function, x)

# Print the result
print(antiderivative)