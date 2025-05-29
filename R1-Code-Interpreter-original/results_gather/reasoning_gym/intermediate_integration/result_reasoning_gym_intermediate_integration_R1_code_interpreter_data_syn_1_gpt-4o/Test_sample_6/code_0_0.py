from sympy import symbols, exp, integrate

# Define the variable and the function
x = symbols('x')
f = 6 * exp(6*x - 6)

# Perform the integration
integral_result = integrate(f, x)

# Print the result
print(integral_result)