import numpy as np

# Step 1: Define the indices for g(r) and S(k) based on the analysis.
# Order: {SS, SR, R, HS, TW}
g_r_indices = [1, 3, 9, 7, 5]
s_k_indices = [4, 0, 2, 8, 6] # Using 0 for the missing S(k) plot for SR

# Step 2: Calculate R_max for the unique system (SR, plot 3).
# Values are estimated from the graph of g(r) for the SR system (Plot 3).
# Based on visual inspection, we assume the calculation should start from r=3/2,
# as g(r) is zero for r<1, making R_g(1/2) undefined.

print("Calculating R_max = max(g(r+1)/g(r)) for the unique system (SR).")
print("Reading values from Plot 3:")
g_1_5 = 0.6 # g(r=1.5)
g_2_5 = 0.8 # g(r=2.5)
g_3_5 = 0.9 # g(r=3.5)

print(f"Value g(r=1.5) = {g_1_5}")
print(f"Value g(r=2.5) = {g_2_5}")
print(f"Value g(r=3.5) = {g_3_5}")

# Calculate the ratios for the first few half-integer r values (starting from r=3/2)
R_g_at_3_div_2 = g_2_5 / g_1_5
R_g_at_5_div_2 = g_3_5 / g_2_5

print("\nCalculating ratios:")
print(f"R_g(3/2) = g(2.5) / g(1.5) = {g_2_5} / {g_1_5} = {R_g_at_3_div_2}")
print(f"R_g(5/2) = g(3.5) / g(2.5) = {g_3_5} / {g_2_5} = {R_g_at_5_div_2}")

# Determine the maximum value
R_max = max(R_g_at_3_div_2, R_g_at_5_div_2)

# Step 3: Combine all values into the final sequence.
final_sequence = g_r_indices + s_k_indices + [R_max]

# Step 4: Print the final answer.
print("\nThe final sequence of 11 values is:")
# Manually formatting to match the desired output style {v1, v2, ...}
sequence_str = "{" + ", ".join(map(str, final_sequence)) + "}"
print(sequence_str)
