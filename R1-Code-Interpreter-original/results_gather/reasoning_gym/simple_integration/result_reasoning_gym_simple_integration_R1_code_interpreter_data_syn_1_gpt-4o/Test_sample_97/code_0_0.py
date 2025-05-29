from sympy import symbols, integrate

# Define the variable
x = symbols('x')

# Define the function
function = -80*x**9 + 25*x**4 + 24*x**3/5

# Calculate the antiderivative
antiderivative = integrate(function, x)

# Print the result
print(antiderivative)