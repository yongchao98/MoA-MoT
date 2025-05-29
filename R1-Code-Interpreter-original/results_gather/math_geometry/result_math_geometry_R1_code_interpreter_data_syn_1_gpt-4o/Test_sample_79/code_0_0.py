import sympy as sp

# Define the variables
r, h = sp.symbols('r h')

# Original volume
V_original = sp.pi * r**2 * h

# New dimensions
r_new = 3 * r
h_new = 2 * h

# New volume
V_new = sp.pi * r_new**2 * h_new

# Calculate the factor X
X = V_new / V_original

# Simplify the expression
X_simplified = sp.simplify(X)

print(X_simplified)