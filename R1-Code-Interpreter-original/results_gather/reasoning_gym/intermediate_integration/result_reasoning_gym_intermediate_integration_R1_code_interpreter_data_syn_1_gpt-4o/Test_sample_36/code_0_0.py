from sympy import symbols, exp, integrate

# Define the variable and the function
x = symbols('x')
f = 16 * exp(-8*x - 6)

# Perform the integration
integral_result = integrate(f, x)

# Print the result
print(integral_result)