from sympy import symbols, integrate, sqrt

# Define the variable and the function
x = symbols('x')
f = 5 * sqrt(-5*x - 4)

# Perform the integration
integral = integrate(f, x)

# Print the result
print(integral)