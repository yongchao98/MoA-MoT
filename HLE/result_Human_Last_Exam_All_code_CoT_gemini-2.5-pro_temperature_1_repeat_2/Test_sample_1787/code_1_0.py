import math

# Constants from the problem statement
c1 = 0.8  # Angular size ratio from the space station
c2 = 0.2  # Angular size ratio from Planet 2

# From the derivation, we assume the ratio of the radii of Planet 2 and the brown dwarf (A)
# is the ratio of the constants c2 and c1.
# A = R2 / Rs
A = c2 / c1
print(f"The derived ratio of the radius of Planet 2 to the brown dwarf's radius (R2/Rs) is: {A}")

# The brightness drop is proportional to the ratio of the disks' areas.
# area_ratio = (pi * R2^2) / (pi * Rs^2) = (R2/Rs)^2
area_ratio = A**2
print(f"The fraction of the brown dwarf's area blocked by Planet 2 is: {area_ratio}")

# The flux during the transit is the initial flux minus the blocked part.
flux_ratio = 1 - area_ratio
print(f"The ratio of flux during transit to the flux before transit is: {flux_ratio}")

# The change in magnitude is calculated using the standard formula.
# delta_m = -2.5 * log10(flux_ratio)
delta_m = -2.5 * math.log10(flux_ratio)

# The final equation is delta_m = -2.5 * log10(1 - (c2/c1)^2)
print(f"The final calculation is: delta_m = -2.5 * log10(1 - ({c2}/{c1})^2)")
print(f"delta_m = -2.5 * log10({flux_ratio})")
print(f"delta_m = {delta_m}")

print(f"\nThe brightness drop of the brown dwarf is {delta_m:.3f} bolometric magnitudes.")

# Final answer format
final_answer = round(delta_m, 3)