import math

# Step 1-4: The ratio of the radii of Planet 2 and the Brown Dwarf.
# As derived from the problem's conditions, this ratio is 0.25.
r2_rs_ratio = 0.25

# Step 5: Calculate the brightness drop.
# The formula for magnitude drop is Δm = -2.5 * log10(1 - (R_p/R_s)^2)
# First, calculate the ratio of the areas, which is (R_p/R_s)^2.
area_ratio = r2_rs_ratio ** 2

# Now, calculate the remaining flux ratio (1 - area_ratio).
flux_ratio = 1 - area_ratio

# Finally, calculate the magnitude drop.
delta_magnitude = -2.5 * math.log10(flux_ratio)

# Print the final equation with the numbers used.
print(f"The calculation for the magnitude drop is:")
print(f"Δm = -2.5 * log10(1 - ({r2_rs_ratio})^2)")
print(f"Δm = -2.5 * log10(1 - {area_ratio})")
print(f"Δm = -2.5 * log10({flux_ratio})")
print(f"Δm = {delta_magnitude:.3f}")

# Final Answer
print(f"\nThe brightness drop is {delta_magnitude:.3f} bolometric magnitudes.")
print(f'<<<__{delta_magnitude:.3f}__>>>')
