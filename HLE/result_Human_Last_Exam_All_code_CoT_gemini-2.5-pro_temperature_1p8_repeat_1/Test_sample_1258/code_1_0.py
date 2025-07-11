import math

# Define the variables symbolically for printing the formula
g = "g"  # Length-to-diameter ratio
k = "k"  # Modulus of the elliptic integrals
E_k = "E(k)"  # Complete elliptic integral of the second kind
F_k = "F(k)"  # Complete elliptic integral of the first kind (as per user's notation)

# The analytical expression for the fluxmetric demagnetizing factor (N_f)
# This formula is derived from the standard expression:
# N_f = (1/g**2) * ( (1/k) * E(k) - (1 - k**2) * F(k) )
# We substitute k = 1 / sqrt(1 + g**2) to express it fully in terms of g.

# Construct the components of the formula as strings
term1_coeff = f"sqrt(1 + {g}**2)"
term2_coeff = f"({g}**2 / (1 + {g}**2))"
expression_in_parentheses = f"{term1_coeff} * {E_k} - {term2_coeff} * {F_k}"
full_formula = f"N_f = (1 / {g}**2) * ({expression_in_parentheses})"

# Print the final result and the definition of terms
print("The analytical expression for the fluxmetric demagnetizing factor (N_f) is:")
print(full_formula)
print("\nWhere:")
print(f"  {g} = the length-to-diameter ratio of the cylinder.")
print(f"  {F_k} = the complete elliptic integral of the first kind with modulus {k}.")
print(f"  {E_k} = the complete elliptic integral of the second kind with modulus {k}.")
print(f"  The modulus {k} is defined as: k = 1 / sqrt(1 + {g}**2)")
