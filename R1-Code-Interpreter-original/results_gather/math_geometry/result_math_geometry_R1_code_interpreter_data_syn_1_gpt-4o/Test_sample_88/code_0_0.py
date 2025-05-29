import sympy as sp

# Define the expression
x = 6 * (sp.sqrt(3) + 1)
AB = 2 * x

# Simplify the expression
AB_simplified = sp.simplify(AB)

# Print the result
print(AB_simplified)