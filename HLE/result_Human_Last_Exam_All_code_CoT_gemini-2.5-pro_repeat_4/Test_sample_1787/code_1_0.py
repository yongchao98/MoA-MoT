import math

# Step 1 & 2: Define constants and derive the relationship between the ratios.
# From the observation at the space station:
# R1/r1 = 0.8 * R2/r2
# From the observation on Planet 2:
# R1/(r2 - r1) = 0.2 * RBD/r2
#
# Let x = r1/r2. Combining these two equations gives:
# R2/RBD = (0.2 / 0.8) * (r2 - r1) / r1 = 0.25 * (1/x - 1) = (1-x)/(4x)

k1 = 0.8  # Angular size ratio from Station
k2 = 0.2  # Angular size ratio from Planet 2

# Step 3: Apply the hidden assumption to find the ratio of orbital radii.
# The assumption is that the radius of the circular orbit (r2) is equal to the
# radius of curvature of the parabolic orbit at its pericenter (r1).
# For a parabolic orbit, the radius of curvature at pericenter is 2 * r1.
# So, r2 = 2 * r1.
# This gives the ratio x = r1/r2.
x = 1.0 / 2.0

print(f"Based on the physical assumption connecting the orbits, the ratio of orbital radii r1/r2 is {x:.1f}.\n")

# Step 4: Solve for the ratio of Planet 2's radius to the Brown Dwarf's radius.
# We use the formula derived from the observational data: R2/RBD = (1-x) / ( (k1/k2) * x )
ratio_k1_k2 = k1 / k2
R2_over_RBD = (1 - x) / (ratio_k1_k2 * x)

print("The final equation for the brightness drop is: Delta_m = -2.5 * log10(1 - (R2/RBD)^2)")
print(f"We need to calculate the value of R2/RBD:")
print(f"R2/RBD = (1 - {x:.1f}) / (({k1}/{k2}) * {x:.1f}) = {1-x} / ({ratio_k1_k2} * {x}) = {R2_over_RBD}\n")

# Step 5: Calculate the brightness drop in magnitudes.
R2_over_RBD_sq = R2_over_RBD**2
flux_ratio = 1 - R2_over_RBD_sq
delta_m = -2.5 * math.log10(flux_ratio)

print(f"Now we calculate the terms inside the magnitude formula:")
print(f"(R2/RBD)^2 = {R2_over_RBD:.4f}^2 = {R2_over_RBD_sq:.4f}")
print(f"1 - (R2/RBD)^2 = 1 - {R2_over_RBD_sq:.4f} = {flux_ratio:.4f}")
print(f"Delta_m = -2.5 * log10({flux_ratio:.4f})")

print(f"\nThe brightness drop of the brown dwarf is {delta_m:.3f} magnitudes.")