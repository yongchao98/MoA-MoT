import math

# Step 1: Determine the ratio of orbital radii, x = r_P1 / r_P2.
# This is derived from the condition that the specific angular momenta are equal.
# l1 = sqrt(2*G*M*r_P1), l2 = sqrt(G*M*r_P2).
# Setting l1 = l2 gives 2*r_P1 = r_P2.
x = 1.0 / 2.0

# Step 2: Determine the ratio of Planet 2's radius to the Brown Dwarf's radius, Z = R_P2 / R_BD.
# This is derived from the two observational equations, leading to the relation:
# 4 * Z * x = 1 - x
# We solve for Z.
Z = (1 - x) / (4 * x)

# Step 3: Calculate the brightness drop in magnitudes, delta_m.
# The formula for brightness drop during a transit is delta_m = -2.5 * log10(1 - (Area_P2 / Area_BD))
# which is delta_m = -2.5 * log10(1 - (R_P2 / R_BD)^2) = -2.5 * log10(1 - Z^2).
delta_m = -2.5 * math.log10(1 - Z**2)

# Output the final equation with the calculated values
print(f"The ratio of the radius of Planet 2 to the radius of the brown dwarf is Z = {Z:.2f}.")
print(f"The brightness drop is calculated using the formula: delta_m = -2.5 * log10(1 - Z^2)")
print(f"delta_m = -2.5 * log10(1 - {Z:.2f}^2)")
print(f"The calculated brightness drop is {delta_m:.3f} magnitudes.")

# Final answer in the required format
print(f"\n<<<The brightness drop is {delta_m:.3f} magnitudes.>>>")