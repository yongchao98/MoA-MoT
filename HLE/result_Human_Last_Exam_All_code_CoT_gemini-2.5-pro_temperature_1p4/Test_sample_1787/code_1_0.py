import math

# Step 1: Define the relationships from the problem statement.
# R_P1 / R_P2 = 0.8 * (r_p1 / r_p2)  (from observation 1)
# R_P1 / R_BD = 0.2 * (1 - r_p1 / r_p2) (from observation 2)
# Let x = r_p1 / r_p2.
# Dividing the two equations or substituting allows us to relate the radius ratios to x.
# (R_P1 / R_BD) / (R_P1 / R_P2) = R_P2 / R_BD
# R_P2 / R_BD = (0.2 * (1 - x)) / (0.8 * x)
# R_P2 / R_BD = 0.25 * (1 - x) / x

# Step 2: Use orbital mechanics to find the value of x = r_p1 / r_p2.
# The problem states P1 is on a parabolic orbit and P2 is on a circular orbit.
# The key physical condition we derive is that their specific angular momentums are equal.
# L1/m1 = sqrt(2*G*M*r_p1) for a parabolic orbit at pericenter.
# L2/m2 = sqrt(G*M*r_p2) for a circular orbit.
# Setting them equal: sqrt(2*G*M*r_p1) = sqrt(G*M*r_p2)
# Squaring both sides: 2*G*M*r_p1 = G*M*r_p2
# 2 * r_p1 = r_p2
# Therefore, x = r_p1 / r_p2 = 1/2
x = 0.5

# Step 3: Substitute x back into the equation for the radius ratio.
r_ratio_p2_bd = 0.25 * (1 - x) / x
print(f"The ratio of the radius of Planet 2 to the Brown Dwarf's radius (R_P2/R_BD) is: {r_ratio_p2_bd}")

# Step 4: Calculate the brightness drop during a transit of Planet 2.
# The ratio of the obscured area to the total area is k = (R_P2 / R_BD)^2
k = r_ratio_p2_bd**2
print(f"The fraction of the Brown Dwarf's area covered by Planet 2 is: {k}")

# The flux ratio during transit is (1 - k)
flux_ratio = 1 - k
print(f"The ratio of flux during transit to normal flux is: {flux_ratio}")

# The brightness drop in magnitudes is delta_m = -2.5 * log10(flux_ratio)
delta_m = -2.5 * math.log10(flux_ratio)

# Final formatted output as requested
print("\nFinal Calculation:")
print(f"Δm = -2.5 * log10(1 - ({r_ratio_p2_bd})^2)")
print(f"Δm = -2.5 * log10({flux_ratio})")
print(f"Δm = {delta_m:.3f} magnitudes")

# The final answer in the required format
final_answer = round(delta_m, 3)
print(f"<<<{final_answer}>>>")