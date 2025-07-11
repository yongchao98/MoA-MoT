# Step 1: Define the parameters based on the problem description.
num_members = 791
num_sections = 61
rows_per_section = num_members // num_sections # This will be 13

initial_radius = 3.0  # r_1 in meters
row_depth = 1.5       # Δr in meters
seated_height = 1.0   # meters
standing_height = 1.5 # meters

# Step 2: Calculate the radii of the first and last rows.
# The radius of row 'i' is given by r_i = initial_radius + (i - 1) * row_depth
first_row_index = 1
last_row_index = rows_per_section

r_1 = initial_radius + (first_row_index - 1) * row_depth
r_13 = initial_radius + (last_row_index - 1) * row_depth

# Step 3: Set up the equation based on the "equality of status" constraint.
# The height of the speaker's head in row 1 is H_speaker = r_1^2/K + standing_height
# The height of the observer's eye in row 13 is H_observer = r_13^2/K + seated_height
# We set H_speaker = H_observer
# r_1^2/K + 1.5 = r_13^2/K + 1.0
# 1.5 - 1.0 = r_13^2/K - r_1^2/K
# 0.5 = (r_13^2 - r_1^2) / K
# K = (r_13^2 - r_1^2) / 0.5
# K = 2 * (r_13^2 - r_1^2)

# Step 4: Calculate K.
K = 2 * (r_13**2 - r_1**2)

# Step 5: Print the final equation with all values and the result.
# We present the equation in the rearranged form K = 2 * (r_13^2 - r_1^2)
print(f"The number of rows per section is {rows_per_section}.")
print(f"The radius of the first row (r_1) is {r_1} m.")
print(f"The radius of the last row (r_13) is {r_13} m.")
print("\nTo satisfy the norm that the speaker in the first row and the observer in the last row have equal eye-level status, we solve the equation:")
print(f"H_speaker = H_observer")
print(f"r_1^2/K + standing_height = r_13^2/K + seated_height")
print(f"{r_1}^2/K + {standing_height} = {r_13}^2/K + {seated_height}")
print("Rearranging for K, we get:")
print(f"K = 2 * (r_13^2 - r_1^2)")
print(f"K = 2 * ({r_13}^2 - {r_1}^2)")
print(f"K = 2 * ({r_13**2} - {r_1**2})")
print(f"K = 2 * ({r_13**2 - r_1**2})")
print(f"K = {int(K)}")

print("\nThe maximum value K can take is an integer.")
print(f"<<<{int(K)}>>>")