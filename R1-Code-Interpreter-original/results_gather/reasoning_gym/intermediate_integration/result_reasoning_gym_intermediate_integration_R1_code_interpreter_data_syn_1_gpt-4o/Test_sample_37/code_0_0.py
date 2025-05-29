from sympy import symbols, exp, diff

# Define the variable and function
x = symbols('x')
f = -3 * exp(10 * x - 6)

# Differentiate the function
f_prime = diff(f, x)

# Print the derivative
print(f_prime)