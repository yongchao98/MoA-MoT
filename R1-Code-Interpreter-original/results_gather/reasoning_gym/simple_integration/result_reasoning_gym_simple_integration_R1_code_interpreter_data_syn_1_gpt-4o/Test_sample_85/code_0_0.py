from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = -30*x**9 + 28*x**6 - x**2/3

# Calculate the antiderivative
antiderivative = integrate(function, x)

# Print the result
print(antiderivative)