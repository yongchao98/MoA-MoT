import math

# Step 1: Define the ratio of the radius of Planet 2 to the Brown Dwarf.
# Based on the problem's geometric and orbital constraints, this ratio (y = R_P2 / R_BD) is determined to be 0.25.
r_p2_div_r_bd = 0.25

# Step 2: Calculate the ratio of fluxes during the transit.
# The drop in flux is proportional to the ratio of the objects' areas.
# Flux_ratio = 1 - (Area_P2 / Area_BD) = 1 - (R_P2 / R_BD)^2
flux_ratio = 1 - r_p2_div_r_bd**2

# Step 3: Calculate the brightness drop in magnitudes.
# Δm = -2.5 * log10(flux_ratio)
delta_m = -2.5 * math.log10(flux_ratio)

# Step 4: Print the final equation with numerical values and the result.
print("The brightness drop (Δm) is calculated using the formula: Δm = -2.5 * log10(1 - (R_P2/R_BD)^2)")
print(f"Plugging in the determined ratio R_P2/R_BD = {r_p2_div_r_bd}:")
print(f"Δm = -2.5 * log10(1 - {r_p2_div_r_bd}^2)")
print(f"Δm = -2.5 * log10(1 - {r_p2_div_r_bd**2})")
print(f"Δm = -2.5 * log10({flux_ratio})")
print(f"The brightness drop is {delta_m:.3f} magnitudes.")

# The final answer in the required format
# print(f'<<<{delta_m:.3f}>>>')