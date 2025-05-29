from sympy import symbols, exp, integrate

# Define the variable and the function
x = symbols('x')
f = (1 - 12*x) * exp(6*x**2 - x - 6)

# Perform the integration
integral = integrate(f, x)

# Output the result
print(integral)