from sympy import symbols, integrate, sin, cos

# Define the variable and the function
x = symbols('x')
u = 3*x + 4
integrand = 3*sin(u)*cos(u)**4

# Perform the integration
integral = integrate(integrand, x)

# Output the result
print(integral)