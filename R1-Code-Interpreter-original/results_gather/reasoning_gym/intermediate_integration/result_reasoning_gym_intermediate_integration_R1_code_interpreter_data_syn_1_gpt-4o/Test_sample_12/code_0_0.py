from sympy import symbols, sin, cos, integrate

# Define the variable and function
x = symbols('x')
u = 4*x + 1
integrand = -4 * sin(u)**2 * cos(u)

# Perform the integration
integral = integrate(integrand, x)

# Output the result
print(integral)