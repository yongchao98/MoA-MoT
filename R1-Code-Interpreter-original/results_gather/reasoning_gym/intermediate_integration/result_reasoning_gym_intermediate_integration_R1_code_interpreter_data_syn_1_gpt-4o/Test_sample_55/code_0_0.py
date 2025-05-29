from sympy import symbols, integrate, sqrt

# Define the variable and the function
x = symbols('x')
f = -15 * sqrt(5*x - 9)

# Perform the integration
integral = integrate(f, x)

# Print the result
print(integral)