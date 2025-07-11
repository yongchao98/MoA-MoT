import math

# Step 1 & 2: Define relationships from the problem statement.
# We derived the ratio of the radii (k = R_2/R_BD) as a function of
# the ratio of the orbital distances (y = d_BD_P2/d_BD_P1).
# k = 0.25 * (y - 1)

# Step 3: Use the hidden constraint from celestial mechanics for this specific alignment.
# The distance ratio for a circular orbit and the pericenter of a parabolic orbit
# in this configuration is 2.
y = 2.0
print(f"The ratio of the orbital distances, d_BD_P2 / d_BD_P1, is assumed to be {y}")

# Step 4: Calculate the ratio of the radii of Planet 2 and the Brown Dwarf.
k = 0.25 * (y - 1)
print(f"The ratio of the radii, R_Planet2 / R_BrownDwarf, is {k}")

# Step 5: Calculate the brightness drop in magnitudes.
# The formula is delta_m = -2.5 * log10(1 - k^2)
k_squared = k**2
flux_ratio = 1 - k_squared
delta_m = -2.5 * math.log10(flux_ratio)

# Output the components of the final equation as requested.
print("\nThe final equation for the brightness drop (delta_m) is:")
print(f"delta_m = -2.5 * log10(1 - (R_Planet2 / R_BrownDwarf)^2)")
print(f"delta_m = -2.5 * log10(1 - {k}^2)")
print(f"delta_m = -2.5 * log10(1 - {k_squared})")
print(f"delta_m = -2.5 * log10({flux_ratio})")
print("\nFinal Answer:")
print(f"The brightness drop is {delta_m:.3f} bolometric magnitudes.")