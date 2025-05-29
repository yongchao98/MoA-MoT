from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
f = -(3*x + 5)**4

# Compute the indefinite integral
integral = integrate(f, x)

# Print the result
print(integral)