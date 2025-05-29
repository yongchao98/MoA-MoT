from sympy import symbols, expand

# Define the variable
y = symbols('y')

# Define the polynomials
poly1 = 2*y**3 + 6*y
poly2 = 23*y**3 + 36*y
poly3 = -19*y**3 - 26*y**2 + 70*y

# Perform the multiplications
result = expand(poly1 * poly2 * poly3)

# Print the simplified result
print(result)