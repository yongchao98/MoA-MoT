import math

# Step 5 & 6: Based on the derived relation and the assumption r2/r1 = 2.
# 1 - r1/r2 = 4 * k * (r1/r2), where k = R2/R_BD
# Let's solve for k
r1_div_r2 = 0.5
# 1 - 0.5 = 4 * k * 0.5
# 0.5 = 2 * k
k = 0.5 / 2
R2_div_R_BD = k

# Step 7: Calculate the flux drop and convert to magnitudes.
# The drop in flux is proportional to the ratio of the areas.
flux_drop_fraction = R2_div_R_BD**2
flux_ratio = 1 - flux_drop_fraction

# The change in magnitude is delta_m = -2.5 * log10(F_new / F_initial)
# The brightness drop is the absolute value of this change.
delta_mag = -2.5 * math.log10(flux_ratio)

# Output the results for the final calculation
print(f"The ratio of the radius of Planet 2 to the Brown Dwarf is R2/R_BD = {R2_div_R_BD}")
print(f"The fraction of the Brown Dwarf's light blocked by Planet 2 is (R2/R_BD)^2 = {flux_drop_fraction}")
print(f"The remaining flux fraction is 1 - {flux_drop_fraction} = {flux_ratio}")
print(f"The brightness drop in magnitudes is -2.5 * log10({flux_ratio:.4f})")
print(f"Final calculated brightness drop: {delta_mag:.3f} magnitudes")

# Final Answer
print("The final equation is:")
print(f"Brightness Drop = -2.5 * log10(1 - ({R2_div_R_BD})**2) = {delta_mag:.3f} magnitudes")