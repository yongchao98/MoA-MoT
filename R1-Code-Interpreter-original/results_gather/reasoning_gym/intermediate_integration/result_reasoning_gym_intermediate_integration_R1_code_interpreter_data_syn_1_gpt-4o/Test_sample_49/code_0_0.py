from sympy import symbols, integrate, cos, sin

# Define the variable and the function
x = symbols('x')
f = -2*x**2*cos(3*x)

# Compute the indefinite integral
integral = integrate(f, x)

# Print the result
print(integral)