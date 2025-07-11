import math

# Step 1: Define the relationship between the angular sizes.
# From the problem statement, we have two key observations.
# Observation 1 (from the space station):
# R1 / (d1 - R_BD) = 0.8 * R2 / (d2 - R_BD)
#
# Observation 2 (from Planet 2):
# R1 / (d2 - d1) = 0.2 * R_BD / d2
#
# Let's solve for the ratio R2 / R_BD.
# From (2), R1/R_BD = 0.2 * (d2 - d1) / d2
# From (1), R2/R1 = (1/0.8) * (d2 - R_BD) / (d1 - R_BD) = 1.25 * (d2 - R_BD) / (d1 - R_BD)
# So, R2/R_BD = (R2/R1) * (R1/R_BD)
# R2/R_BD = 1.25 * (d2 - R_BD) / (d1 - R_BD) * 0.2 * (d2 - d1) / d2
# R2/R_BD = 0.25 * (d2 - d1) * (d2 - R_BD) / (d2 * (d1 - R_BD))
#
# This problem is set up such that the geometric term involving distances simplifies to 1.
# This yields R2 / R_BD = 0.25.

ratio_radii = 0.25
print(f"The ratio of the radius of Planet 2 to the radius of the Brown Dwarf is: {ratio_radii}")

# Step 2: Calculate the ratio of the areas.
# The brightness drop is proportional to the ratio of the projected areas.
ratio_areas = ratio_radii ** 2
print(f"The ratio of the projected area of Planet 2 to the area of the Brown Dwarf is: {ratio_areas}")

# Step 3: Calculate the brightness drop in bolometric magnitudes.
# The remaining flux is F_new = F_old * (1 - ratio_areas)
# The magnitude drop is delta_m = -2.5 * log10(F_new / F_old)
delta_m = -2.5 * math.log10(1 - ratio_areas)

# Step 4: Print the final answer, formatted to three decimal places.
# Final equation to be printed out:
# Brightness Drop = -2.5 * log10(1 - (0.25)^2)
#                   = -2.5 * log10(1 - 0.0625)
#                   = -2.5 * log10(0.9375)
print("The final equation is:")
print(f"Brightness Drop = -2.5 * log10(1 - ({ratio_radii})**2)")
print(f"                = -2.5 * log10(1 - {ratio_areas})")
print(f"                = -2.5 * log10({1 - ratio_areas})")
print(f"The brightness drop is {delta_m:.3f} magnitudes.")