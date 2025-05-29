from sympy import symbols, integrate, cos, sin

# Define the variable and the function
x = symbols('x')
f = -2*x*cos(x)

# Perform the integration
integral_result = integrate(f, x)

# Print the result
print(integral_result)