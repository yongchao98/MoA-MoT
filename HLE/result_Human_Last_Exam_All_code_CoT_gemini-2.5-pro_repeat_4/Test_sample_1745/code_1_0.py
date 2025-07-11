# The problem requires finding the expression for the Electrical Double-Layer (EDL)
# potential distribution, denoted as ψ(y), in a parallel-plate microchannel.
#
# Based on the linearized Poisson-Boltzmann equation and the specified boundary conditions:
# 1. At the top wall (y = H/2), the potential ψ is 0.
# 2. At the bottom wall (y = -H/2), the potential ψ is z₁(1 + βk).
#
# The final derived expression for ψ(y) is presented below. The equation includes the
# numbers 1 and 2 as part of the terms z₁ and H/2, respectively.

# Define the components of the equation as symbolic strings
psi_of_y = "ψ(y)"
zeta_potential_term = "z₁(1 + βk)"
numerator = "sinh(k(H/2 - y))"
denominator = "sinh(kH)"

# Construct the final equation string
final_equation = f"{psi_of_y} = {zeta_potential_term} * ({numerator} / {denominator})"

# Print the final expression
print("The expression for the Electrical double-layer potential distribution is:")
print(final_equation)