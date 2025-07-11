# Based on the analysis of the plots, the indices for each simulation are determined.

# g(r) plot indices for the sequence {SS, SR, R, HS, TW}
g_ss_idx = 3
g_sr_idx = 7
g_r_idx = 1
g_hs_idx = 4
g_tw_idx = 5

# S(k) plot indices for the sequence {SS, SR, R, HS, TW}
# The Ramp (R) system is identified as the unique system with no S(k) plot.
s_ss_idx = 9
s_sr_idx = 2
s_r_idx = 0  # Not plotted
s_hs_idx = 8
s_tw_idx = 6

# The unique system is Ramp (R), and its g(r) is Plot 1.
# R_max is the maximum of g(r+1)/g(r) for r = {1.5, 2.5, 3.5, ...} (dimensionless r/sigma).
# Values are estimated from Plot 1.
# At r = 2.5, the ratio g(3.5)/g(2.5) gives the maximum.
g_at_2_5 = 0.9
g_at_3_5 = 1.05
r_max = g_at_3_5 / g_at_2_5

# The final answer is a sequence of 11 values.
final_sequence = f"{{{g_ss_idx}, {g_sr_idx}, {g_r_idx}, {g_hs_idx}, {g_tw_idx}, {s_ss_idx}, {s_sr_idx}, {s_r_idx}, {s_hs_idx}, {s_tw_idx}, {r_max}}}"

print(final_sequence)