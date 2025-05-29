from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
f = -2 * (8*x - 7)**4

# Perform the integration
integral_result = integrate(f, x)

# Print the result
print(integral_result)