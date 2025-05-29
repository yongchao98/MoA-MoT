from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
f = 2 * (3*x + 10)**3

# Compute the indefinite integral
integral = integrate(f, x)

# Print the result
print(integral)