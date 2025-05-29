from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
f = 3 * (4*x - 8)**3

# Perform the integration
integral = integrate(f, x)

# Print the result
print(integral)