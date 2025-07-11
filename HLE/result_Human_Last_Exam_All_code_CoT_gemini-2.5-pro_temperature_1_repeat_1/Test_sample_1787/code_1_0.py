import math

# Step 1 & 2: Define ratios from the problem statement
# k1 is the ratio of angular sizes of P1 to P2 as seen from the station.
k1 = 0.8
# k2 is the ratio of angular sizes of P1 to the Brown Dwarf as seen from P2.
k2 = 0.2

# Step 3: Make a simplifying assumption
# The problem is under-determined. We assume a special geometric configuration
# where Planet 1 is equidistant from the Brown Dwarf and Planet 2.
# This means d_BD1 = 0.5 * d_BD2, so the ratio x = d_BD1 / d_BD2 is 0.5.
x = 0.5
print(f"Assuming Planet 1 is halfway between the Brown Dwarf and Planet 2, the ratio of their orbital distances is x = {x}")

# Step 4: Calculate the radius ratio z = R_2 / R_BD
# From the derived relation z = (1 - x) / (4 * x)
z = (1 - x) / (4 * x)
print(f"The calculated ratio of Planet 2's radius to the Brown Dwarf's radius is z = {z}")

# Step 5: Calculate the brightness drop in magnitudes
# The formula for the brightness drop (delta_m) is:
# delta_m = -2.5 * log10(1 - (R_2/R_BD)^2)
# which is -2.5 * log10(1 - z^2)

# Calculate the term inside the logarithm
flux_ratio = 1 - z**2

# Calculate the final magnitude drop
delta_m = -2.5 * math.log10(flux_ratio)

# Print the final equation with all the numbers
print(f"\nThe final equation for the brightness drop is:")
print(f"Δm = -2.5 * log10(1 - ({z})^2)")
print(f"Δm = -2.5 * log10({flux_ratio})")

# Print the result rounded to three decimal places
print(f"\nThe brightness drop of the brown dwarf is {delta_m:.3f} bolometric magnitudes.")

# Final answer in the specified format
print(f"\n<<<{delta_m:.3f}>>>")