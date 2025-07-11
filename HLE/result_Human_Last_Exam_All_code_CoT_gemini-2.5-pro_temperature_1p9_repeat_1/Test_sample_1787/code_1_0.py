import math

# Step 1 & 2: Derive the relationship between radius and distance ratios.
# From the problem statement, we have:
# Obs 1: (R1/d1) = 0.8 * (R2/d2)
# Obs 2: R1/(d2-d1) = 0.2 * (R_BD/d2)
#
# Combining these, we get:
# 4 * (R2 / R_BD) = (d2 / d1) - 1
# This is our first key equation relating the radius ratio to the distance ratio.

# Step 3: Use orbital mechanics to find the distance ratio.
# The specific angular momentum (h) of a body in orbit around a central mass M is h = r * v.
# For Planet 1 at the pericenter (d1) of a parabolic orbit: v1 = sqrt(2*G*M/d1)
# So, h1 = d1 * sqrt(2*G*M/d1) = sqrt(2*G*M*d1)
# For Planet 2 in a circular orbit (d2): v2 = sqrt(G*M/d2)
# So, h2 = d2 * sqrt(G*M/d2) = sqrt(G*M*d2)
#
# A key physical assumption for this configuration is that their specific angular momenta are equal (h1 = h2).
# sqrt(2*G*M*d1) = sqrt(G*M*d2)
# Squaring both sides gives: 2*G*M*d1 = G*M*d2
# This simplifies to d2 = 2*d1, so the distance ratio d2/d1 = 2.
d2_div_d1 = 2

# Step 4: Solve for the radius ratio R2/R_BD.
# Substitute d2/d1 = 2 into our first key equation:
# 4 * (R2 / R_BD) = 2 - 1 = 1
# So, R2 / R_BD = 1/4
R2_div_R_BD = 0.25

# Step 5: Calculate the brightness drop factor.
# This is the ratio of the areas, (R2/R_BD)^2.
brightness_drop_factor = R2_div_R_BD ** 2

# Step 6: Calculate the change in magnitude.
# The formula is Δm = -2.5 * log10(1 - brightness_drop_factor).
factor_const = -2.5
base_val = 1
flux_ratio = base_val - brightness_drop_factor
delta_m = factor_const * math.log10(flux_ratio)

# Output the components of the final equation and the result.
print("The final equation for the change in magnitude (Δm) is:")
print("Δm = C * log10(B - F)")
print(f"Where C = {factor_const}, B = {base_val}, and F (the brightness drop factor (R_planet/R_star)^2) = {brightness_drop_factor}")
print("\nCalculating the result:")
print(f"Δm = {factor_const} * log10({base_val} - {brightness_drop_factor})")
print(f"Δm = {factor_const} * log10({flux_ratio})")
print(f"The brightness drop in bolometric magnitudes is: {delta_m:.3f}")

# Final Answer in the required format
final_answer = round(delta_m, 3)
print(f"\n<<<{final_answer}>>>")