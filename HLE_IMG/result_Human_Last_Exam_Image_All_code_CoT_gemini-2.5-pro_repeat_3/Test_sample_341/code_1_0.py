import math

# Step 1: Define the plot indices based on the analysis.
# g(r) for {SS, SR, R, HS, TW}
g_indices = [1, 3, 5, 7, 9]
# S(k) for {SS, SR, R, HS, TW}
s_indices = [8, 6, 0, 2, 4]

# Step 2: Define the calculated value for R_max.
# The unique system is R (Ramp), its g(r) is Plot 5.
# R_max = max_{r=3/2, 5/2, ...} g(r+1)/g(r)
# From Plot 5, g(1.5)=1.0, g(2.5)=0.8, g(3.5)=0.9
# R_g(1.5) = g(2.5)/g(1.5) = 0.8/1.0 = 0.8
# R_g(2.5) = g(3.5)/g(2.5) = 0.9/0.8 = 1.125
R_max = 1.125

# Step 3: Combine all values into a single list.
final_values = g_indices + s_indices + [R_max]

# Step 4: Print the final result in the specified format.
# "output each number in the final equation!"
print(f"{{{final_values[0]}, {final_values[1]}, {final_values[2]}, {final_values[3]}, {final_values[4]}, {final_values[5]}, {final_values[6]}, {final_values[7]}, {final_values[8]}, {final_values[9]}, {final_values[10]}}}")
