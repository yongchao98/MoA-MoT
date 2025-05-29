from sympy import symbols, integrate, sqrt

# Define the variable and the function
x = symbols('x')
f = 12*sqrt(3 - 4*x)

# Perform the integration
integral = integrate(f, x)

# Output the result
print(integral)