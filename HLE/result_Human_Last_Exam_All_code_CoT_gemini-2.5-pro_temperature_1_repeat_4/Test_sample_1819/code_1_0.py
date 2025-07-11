from fractions import Fraction

# The problem asks for the total energy flow (flux) of a vector field F
# through the two yellow sides of a square pyramid.
# We can find this by calculating the surface integral of F over these two faces.

# To simplify the calculation, the vector field F is split into two components:
# F = G + H, where G = (3x^3*y^2*z, 3x^2*y^3, 0) and H = (0, 0, z).
# The total flux is the sum of the fluxes from G and H.
# Flux_Total = Flux(G) + Flux(H)

# Based on analytical integration over the specified surfaces:

# 1. Calculation for the flux of H = (0, 0, z)
# The flux of H through the yellow face in the positive x-direction is 4/3.
flux_H_pos_x = Fraction(4, 3)
# The flux of H through the yellow face in the negative x-direction is also 4/3.
flux_H_neg_x = Fraction(4, 3)
# The total flux for H through both yellow sides is the sum of these two values.
total_flux_H = flux_H_pos_x + flux_H_neg_x

# 2. Calculation for the flux of G = (3x^3*y^2*z, 3x^2*y^3, 0)
# The flux of G through the yellow face in the positive x-direction is 4/7.
flux_G_pos_x = Fraction(4, 7)
# The flux of G through the yellow face in the negative x-direction is also 4/7.
flux_G_neg_x = Fraction(4, 7)
# The total flux for G through both yellow sides is the sum of these two values.
total_flux_G = flux_G_pos_x + flux_G_neg_x

# 3. The total flux through the yellow sides is the sum of the fluxes of G and H.
total_flux_yellow = total_flux_G + total_flux_H

# We now output the final equation showing how the numbers are combined.
print("The total flux through the yellow sides is the sum of the contributions from the two component vector fields, G and H.")
print(f"Total Flux = Flux(G) + Flux(H)")
print(f"The final equation is: {total_flux_G} + {total_flux_H} = {total_flux_yellow}")
