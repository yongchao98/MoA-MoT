import sympy as sp

# Define the variable and the function
x = sp.symbols('x')
f = -2*x**3*sp.sin(x)

# Compute the indefinite integral
integral = sp.integrate(f, x)

# Print the result
print(integral)