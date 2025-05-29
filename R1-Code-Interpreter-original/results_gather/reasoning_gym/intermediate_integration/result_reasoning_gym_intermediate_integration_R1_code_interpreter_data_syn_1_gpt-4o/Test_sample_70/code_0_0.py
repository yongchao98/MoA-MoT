from sympy import symbols, integrate, sqrt

# Define the variable and the function
x = symbols('x')
f = -18 * sqrt(9*x + 7)

# Perform the integration
integral = integrate(f, x)

# Print the result
print(integral)