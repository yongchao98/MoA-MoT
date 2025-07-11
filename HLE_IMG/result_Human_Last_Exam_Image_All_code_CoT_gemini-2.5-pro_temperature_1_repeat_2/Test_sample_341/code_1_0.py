import collections

# Step 1: Define the plot indices based on the analysis.
# The order of systems is {SS, SR, R, HS, TW}.
g_r_indices = collections.OrderedDict([
    ('SS', 1),
    ('SR', 3),
    ('R', 7),
    ('HS', 9),
    ('TW', 5)
])

s_k_indices = collections.OrderedDict([
    ('SS', 2),
    ('SR', 6),
    ('R', 0),  # Not plotted
    ('HS', 4),
    ('TW', 8)
])

# Step 2: Calculate R_max.
# The unique system is R, corresponding to g(r) plot 7.
# R_g(r) = g(r+1)/g(r). We need to evaluate this for r = {1/2, 3/2, 5/2, ...}.
# For r=1/2, g(1/2)=0, so the ratio is undefined. We start with r=3/2.
# From visual inspection of plot 7:
# g(r=1.5) is at the first minimum, approx 0.7.
# g(r=2.5) is at the second minimum, approx 0.9.
# The oscillations are damping, so the maximum ratio is expected for the first valid r.
g_at_1_5 = 0.7
g_at_2_5 = 0.9
R_max = g_at_2_5 / g_at_1_5

# Step 3: Assemble the final list of 11 values.
final_values = list(g_r_indices.values()) + list(s_k_indices.values()) + [R_max]

# Step 4: Print the final answer in the required format.
# The problem asks to output each number in the final result.
# We will format the R_max value to a reasonable number of decimal places.
output_string = "{"
output_string += str(g_r_indices['SS']) + ", "
output_string += str(g_r_indices['SR']) + ", "
output_string += str(g_r_indices['R']) + ", "
output_string += str(g_r_indices['HS']) + ", "
output_string += str(g_r_indices['TW']) + ", "
output_string += str(s_k_indices['SS']) + ", "
output_string += str(s_k_indices['SR']) + ", "
output_string += str(s_k_indices['R']) + ", "
output_string += str(s_k_indices['HS']) + ", "
output_string += str(s_k_indices['TW']) + ", "
output_string += f"{R_max}"
output_string += "}"

print(output_string)