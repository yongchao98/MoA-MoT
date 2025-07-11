import numpy as np

# Step 1: Define the plot indices based on the analysis.
# Order of systems: SS, SR, R, HS, TW
g_indices = [1, 9, 3, 7, 5]
s_indices = [6, 4, 0, 8, 2]

# Step 2: Calculate the 11th value, R_max, for the unique system (Ramp, plot 3).
# The calculation is R_max = g(r+1)/g(r) for r=1.5, as this gives the maximum value.
# The values are read from the linear segments in plot 3, assuming the blue dotted line is at y=1.
# On the segment from r=1 to r=2, the function goes from ~1.8 to 1.0. At the midpoint r=1.5, g(1.5)=1.4.
# On the segment from r=2 to r=3, the function goes from ~0.6 to 0.2. At the midpoint r=2.5, g(2.5)=0.4.
g_at_1_5 = 1.4
g_at_2_5 = 0.4

# Perform the calculation for R_max.
R_max = g_at_2_5 / g_at_1_5

# As requested, output the numbers in the final equation for R_max.
print("Calculation of the 11th value, R_max:")
print(f"R_max = g(1.5 + 1) / g(1.5)")
print(f"Using values from Plot 3: g(1.5) = {g_at_1_5}, g(2.5) = {g_at_2_5}")
print(f"R_max = {g_at_2_5} / {g_at_1_5} = {R_max}")

# Step 3: Combine all values into the final answer string for the user.
all_values = g_indices + s_indices + [R_max]
answer_string = "{" + ",".join([str(v) for v in all_values]) + "}"

print("\nFinal list of 11 values:")
print(answer_string)