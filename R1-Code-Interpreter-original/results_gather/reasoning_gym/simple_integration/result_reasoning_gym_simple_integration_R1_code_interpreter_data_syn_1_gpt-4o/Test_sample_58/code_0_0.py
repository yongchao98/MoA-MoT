from sympy import symbols, integrate

# Define the variable
x = symbols('x')

# Define the function
function = -8*x**7 - 7*x**6 + 10*x**5/3

# Calculate the antiderivative
antiderivative = integrate(function, x)

# Print the result
print(antiderivative)