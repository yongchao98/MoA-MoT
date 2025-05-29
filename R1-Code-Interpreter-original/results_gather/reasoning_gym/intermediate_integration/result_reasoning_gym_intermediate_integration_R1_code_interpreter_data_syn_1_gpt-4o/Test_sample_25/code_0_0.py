from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
f = -2 * (4*x - 9)**3

# Perform the integration
integral = integrate(f, x)

# Output the result
print(integral)