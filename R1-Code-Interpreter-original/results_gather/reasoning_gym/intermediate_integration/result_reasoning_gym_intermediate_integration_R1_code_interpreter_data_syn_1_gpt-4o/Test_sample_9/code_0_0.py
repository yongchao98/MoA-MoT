from sympy import symbols, integrate, sin, cos

# Define the variable and function
x = symbols('x')
u = 10*x + 9
integrand = 20*sin(u)**3*cos(u)

# Perform the integration
integral = integrate(integrand, x)

# Print the result
print(integral)