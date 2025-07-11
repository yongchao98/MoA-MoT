import sympy

# Define the symbols for the variables in the problem.
# mu_0 is the permeability of free space.
# h is the separation between wires in a circuit.
# d is the distance between the two circuits.
# R1 is the inner radius of the concentrator shells.
# R2 is the outer radius of the concentrator shells.
mu_0, h, d, R1, R2, pi = sympy.symbols('mu_0 h d R1 R2 pi')

# The mutual inductance per unit length for the bare circuits is M1/L.
# M1_per_L = (mu_0 * h**2) / (2 * pi * d**2)
# The mutual inductance per unit length with the concentrators is M2/L.
# M2_per_L = M1_per_L * (R2/R1)**2
# The change in mutual inductance per unit length is Delta_M/L = M2/L - M1/L.
delta_M_expr = (mu_0 * h**2) / (2 * pi * d**2) * ((R2/R1)**2 - 1)

# To present the formula in a readable string format
# We build the string representing the formula step by step
term1 = "(mu_0 * h**2)"
term2 = "(2 * pi * d**2)"
term3 = "((R2/R1)**2 - 1)"

print("The expression for the change in mutual inductance per unit length (Delta M / L) is:")
print(f"Delta_M_per_L = {term1} / {term2} * {term3}")

# Let's print the symbolic expression object from sympy as well, for clarity
# Note: The prompt asks for the expression with numbers, the above string formatting handles this best.
# The below print is for additional verification and clarity of the object used.
# from sympy import pprint
# print("\nSympy expression object:")
# pprint(delta_M_expr)
