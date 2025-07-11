import math

# Step 1 & 2: Identify the plots based on physical characteristics.
# The order of systems is {SS, SR, R, HS, TW}.

# g(r) indices:
# SS (Plot 1): Repulsive shoulder suppresses g(r) for 1 < r < 1.5.
# SR (Plot 3): Sticky attraction causes a delta-function-like peak at r=1.
# R (Plot 5): Repulsive ramp causes g(r) to increase linearly for 1 < r < 1.5.
# HS (Plot 7): Standard hard-sphere g(r) with a finite contact peak.
# TW (Plot 9): Attractive well enhances g(r) for 1 < r < 1.5.
g_r_indices = [1, 3, 5, 7, 9]

# S(k) indices:
# HS (Plot 6): S(0) = (1-1/3)^2 ~ 0.44, matching the plot.
# SR (Plot 4): Strong attraction leads to very high S(0).
# SS (Plot 2): Strong repulsion leads to very low S(0).
# R (Plot 8): Weaker repulsion leads to a low S(0), but higher than SS.
# TW: No S(k) plot remains. This is the unique system.
s_k_indices = [2, 4, 8, 6, 0]

# Step 3: The unique system is Triangle Well (TW), with g(r) shown in Plot 9.

# Step 4: Calculate R_max for the TW system.
# R_g(r) = g(r+1)/g(r) for r = {1/2, 3/2, 5/2, ...}.
# We start with r=3/2 since g(1/2)=0.
# Values are read from Plot 9 (g(r) for TW), where the x-axis is r/sigma.
# We estimate the values from the grid, where major lines are at intervals of 0.5.

# At r = 1.5: g(1.5) is at the end of the attractive well. g(1.5) ≈ 0.8
g_1_5 = 0.8
# At r = 2.5: g(2.5) is the peak of the first oscillation. g(2.5) ≈ 1.1
g_2_5 = 1.1
# At r = 3.5: g(3.5) is the peak of the second oscillation. g(3.5) ≈ 1.05
g_3_5 = 1.05

# Calculate the first ratio for r = 1.5:
# R_g(1.5) = g(1.5 + 1) / g(1.5) = g(2.5) / g(1.5)
r_g_1_5 = g_2_5 / g_1_5

# Calculate the second ratio for r = 2.5:
# R_g(2.5) = g(2.5 + 1) / g(2.5) = g(3.5) / g(2.5)
r_g_2_5 = g_3_5 / g_2_5

# The oscillations in g(r) are damping towards 1, so the ratio g(r+1)/g(r) will also
# approach 1. The maximum value is the first one calculated.
R_max = r_g_1_5

# Step 5: Assemble the final list of 11 values.
final_values = g_r_indices + s_k_indices + [R_max]

# Print the final result in the specified format.
# The problem asks to output each number in the final equation.
g_ss, g_sr, g_r, g_hs, g_tw = g_r_indices
s_ss, s_sr, s_r, s_hs, s_tw = s_k_indices
print(f"The identified plots for g(r) in order {{SS, SR, R, HS, TW}} are: {{{g_ss}, {g_sr}, {g_r}, {g_hs}, {g_tw}}}")
print(f"The identified plots for S(k) in order {{SS, SR, R, HS, TW}} are: {{{s_ss}, {s_sr}, {s_r}, {s_hs}, {s_tw}}}")
print(f"The unique system is TW (plot {g_tw}).")
print(f"To calculate R_max, we evaluate R_g(r) = g(r+1)/g(r) for r = 1.5, 2.5, ...")
print(f"From plot {g_tw}: g(1.5) = {g_1_5}, g(2.5) = {g_2_5}")
print(f"R_g(1.5) = g(2.5) / g(1.5) = {g_2_5} / {g_1_5} = {r_g_1_5}")
print(f"From plot {g_tw}: g(3.5) = {g_3_5}")
print(f"R_g(2.5) = g(3.5) / g(2.5) = {g_3_5} / {g_2_5} = {r_g_2_5:.3f}")
print(f"The maximum value is R_max = {R_max}")
print("\nThe final sequence is:")
print(f"{{{','.join(map(str, final_values))}}}")