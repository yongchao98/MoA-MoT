# Step 1: Define the identified plot indices for g(r)
# Sequence is {SS, SR, R, HS, TW}
g_r_indices = {
    "SS": 1,
    "SR": 3,
    "R": 5,
    "HS": 7,
    "TW": 9
}
g_r_sequence = [g_r_indices["SS"], g_r_indices["SR"], g_r_indices["R"], g_r_indices["HS"], g_r_indices["TW"]]

# Step 2: Define the identified plot indices for S(k)
# Sequence is {SS, SR, R, HS, TW}. HS is the unique system with no S(k) plot (index 0).
s_k_indices = {
    "SS": 6,
    "SR": 2,
    "R": 8,
    "HS": 0,
    "TW": 4
}
s_k_sequence = [s_k_indices["SS"], s_k_indices["SR"], s_k_indices["R"], s_k_indices["HS"], s_k_indices["TW"]]

# Step 3: Calculate R_max for the unique system (HS, plot 7)
# R_g(r) = g(r+1)/g(r)
# We read values from plot 7 for r = 1.5, 2.5, ...
# For r=1.5:
g_at_1_5 = 1.0
g_at_2_5 = 1.1
R_g_1_5 = g_at_2_5 / g_at_1_5

# For r=2.5:
g_at_3_5 = 1.05
R_g_2_5 = g_at_3_5 / g_at_2_5

# The function g(r) approaches 1, so the ratios R_g(r) also approach 1.
# The maximum value is the first one.
R_max = R_g_1_5

# Step 4: Assemble and print the final result in the required format
final_sequence = g_r_sequence + s_k_sequence + [R_max]

# The problem asks to output each number in the final equation.
# We will format the output string as {val1, val2, ..., val11}.
output_str = "{" + ", ".join(map(str, final_sequence)) + "}"
print(output_str)