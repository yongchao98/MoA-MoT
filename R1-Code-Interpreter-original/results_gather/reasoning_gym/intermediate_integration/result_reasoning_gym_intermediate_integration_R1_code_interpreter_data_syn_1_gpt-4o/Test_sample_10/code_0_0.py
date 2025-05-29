from sympy import symbols, exp, integrate

# Define the variable and function
x = symbols('x')
u = 4*x**2 + 2*x + 10
f = (16*x + 4) * exp(u)

# Perform the integration
integral_result = integrate(f, x)

# Print the result
print(integral_result)