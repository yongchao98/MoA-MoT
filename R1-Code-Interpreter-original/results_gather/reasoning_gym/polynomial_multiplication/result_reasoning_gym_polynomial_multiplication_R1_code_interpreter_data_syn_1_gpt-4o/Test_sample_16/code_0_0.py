from sympy import symbols, expand

# Define the variable
z = symbols('z')

# Define the polynomials
P1 = 51*z**2 - 85
P2 = 72*z**3 + 65*z**2
P3 = -66*z**2 + 17*z + 13

# Multiply the polynomials
result = expand(P1 * P2 * P3)

# Print the expanded and simplified result
print(result)