from sympy import symbols, integrate, sin, cos

# Define the variable and the function
x = symbols('x')
f = -sin(x + 6) * cos(x + 6)**4

# Perform the integration
integral_result = integrate(f, x)

# Print the result
print(integral_result)