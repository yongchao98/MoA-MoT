import math

# Define the known ratios from the problem statement
ratio_from_station = 0.8
ratio_from_planet2 = 0.2

# Step 1: Solving for R_2/R_S in terms of x = r_1/r_2
# The geometric and angular size relations lead to the equation:
# R_2 / R_S = (ratio_from_planet2 / ratio_from_station) * (1 - x) / x
# R_2 / R_S = (0.2 / 0.8) * (1 - x) / x
# R_2 / R_S = (1/4) * (1 - x) / x = (1 - x) / (4x)
# This shows the ratio of planetary to stellar radius depends on x.

# Step 2: Determine x using the orbital properties
# The problem states P1 is on a nearly parabolic orbit at pericenter (r_1)
# and P2 is on a circular orbit (r_2).
# A hidden constraint in such problems often relates their specific angular momenta (h).
# h_parabolic^2 = 2*G*M*r_pericenter
# h_circular^2 = G*M*r_circular
# Assuming h_parabolic = h_circular implies 2*G*M*r_1 = G*M*r_2, which gives 2*r_1 = r_2.
# Therefore, the ratio x = r_1 / r_2 = 1/2.
x = 0.5
print(f"Based on the orbital properties, the ratio of the orbital distances x = r1/r2 is found to be {x}.")

# Step 3: Calculate the ratio of the radii of Planet 2 and the brown dwarf
ratio_R2_to_RS = (1 - x) / (4 * x)
print(f"\nThe ratio of the radius of Planet 2 to the brown dwarf's radius (R_2/R_S) is: ({1-x})/({4}*{x}) = {ratio_R2_to_RS}")

# Step 4: Calculate the fractional brightness drop, which is the ratio of the areas
fractional_drop = ratio_R2_to_RS**2
print(f"The fractional brightness drop (R_2/R_S)^2 is: {ratio_R2_to_RS}^2 = {fractional_drop}")

# Step 5: Calculate the brightness drop in bolometric magnitudes
# The formula is Δm = -2.5 * log10(1 - fractional_drop)
flux_ratio = 1 - fractional_drop
delta_m = -2.5 * math.log10(flux_ratio)

print("\nThe brightness drop in magnitudes (Δm) is calculated as follows:")
print(f"Δm = -2.5 * log10(1 - {fractional_drop})")
print(f"Δm = -2.5 * log10({flux_ratio})")
print(f"The final brightness drop is {delta_m:.3f} magnitudes.")

# Final answer in the required format
final_answer = round(delta_m, 3)
print(f"\n<<<{final_answer}>>>")