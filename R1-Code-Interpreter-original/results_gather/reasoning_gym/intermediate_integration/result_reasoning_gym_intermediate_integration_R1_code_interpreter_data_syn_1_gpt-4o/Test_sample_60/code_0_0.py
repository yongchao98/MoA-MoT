from sympy import symbols, sin, cos, integrate

# Define the variable and function
x = symbols('x')
u = 3*x - 10
integrand = 3*sin(u)**3*cos(u)

# Perform the integration
integral = integrate(integrand, x)

# Output the result
print(integral)