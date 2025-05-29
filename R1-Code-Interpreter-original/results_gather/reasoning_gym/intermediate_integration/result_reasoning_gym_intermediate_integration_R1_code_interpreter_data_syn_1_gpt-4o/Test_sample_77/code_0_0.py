from sympy import symbols, asin, integrate, sqrt

# Define the variable and function
x = symbols('x')
integrand = asin(3*x)

# Perform the integration
integral_result = integrate(integrand, x)

# Output the result
print(integral_result)