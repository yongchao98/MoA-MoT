from sympy import symbols, sin, integrate

# Define the variable and function
x = symbols('x')
u = 4*x - 6
integrand = -12*sin(u)**3*sin(u).diff(x)

# Perform the integration
integral = integrate(integrand, x)

# Output the result
print(integral)