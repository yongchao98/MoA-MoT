from sympy import symbols, integrate, sin, cos

# Define the variable and the function
x = symbols('x')
f = -10*sin(5*x + 3)*cos(5*x + 3)

# Perform the integration
integral_result = integrate(f, x)

# Print the result
print(integral_result)