import math

# Step 1: Define the ratio of radii based on the orbital mechanics insight.
# The ratio of Planet 2's orbital radius to Planet 1's pericenter radius is 2.
r2_div_r1 = 2

# Step 2: Calculate the ratio of the planet's radius to the brown dwarf's radius.
# From the angular size relationships, we derived R_2 / R_BD = 0.25 * (r_2/r_1 - 1).
r_p2_div_r_bd = 0.25 * (r2_div_r1 - 1)

# Step 3: Calculate the flux drop ratio.
# Flux drop is the square of the ratio of the radii.
flux_drop_ratio = r_p2_div_r_bd**2

# Step 4: Calculate the brightness drop in magnitudes.
# The formula is delta_m = -2.5 * log10(1 - flux_drop_ratio).
delta_m = -2.5 * math.log10(1 - flux_drop_ratio)

# Step 5: Print the final equation with the calculated values
# to show the calculation steps as requested.
print("The final equation for the brightness drop in magnitudes (Δm) is:")
print(f"Δm = -2.5 * log10(1 - ({r_p2_div_r_bd:.2f})**2)")
print(f"Δm = -2.5 * log10(1 - {flux_drop_ratio})")
print(f"Δm = -2.5 * log10({1 - flux_drop_ratio})")

# Step 6: Print the final answer, rounded to three decimal places.
print("\nThe calculated brightness drop is:")
print(f"{delta_m:.3f}")

print("\nFinal Answer in desired format:")
print(f"<<<{delta_m:.3f}>>>")