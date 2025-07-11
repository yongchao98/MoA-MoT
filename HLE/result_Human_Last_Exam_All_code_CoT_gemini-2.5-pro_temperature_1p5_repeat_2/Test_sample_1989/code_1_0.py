import sympy

# Define the variables
r, theta, A, B = sympy.symbols('r theta A B')

# The corrector term K(theta)
K = A * (1 - sympy.cos(theta)) + B * sympy.sin(theta)

# The full corrector is r**K
corrector = r**K

# Print the corrector term as a string
# This represents the multiplicative factor to the solution for A=B=0
print(f"{corrector}")