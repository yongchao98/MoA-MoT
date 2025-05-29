from sympy import symbols, integrate, exp

# Define the variable and the function
x = symbols('x')
f = -3*x**3*exp(x)

# Compute the indefinite integral
integral_result = integrate(f, x)

# Print the result
print(integral_result)