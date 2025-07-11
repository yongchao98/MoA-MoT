import math

# Step 1: Define the parameters based on the problem description.
total_members = 791
num_sections = 61
initial_radius = 3.0  # r_1 in meters
row_depth = 1.5       # distance between rows in meters
speaker_height_add = 1.5 # height of a standing person from the floor
seated_height_add = 1.0  # height of a seated person's head/eyes from the floor

# Step 2: Calculate the number of rows.
num_rows = total_members // num_sections

# Step 3 & 4: Formulate the visibility constraint and find the most restrictive case.
# The visibility constraint, when solved for K, gives:
# K < 2 * (r_i - r_1) * (r_j - r_1)
# To satisfy this for all i and j, K must be less than the minimum possible value of the right side.
# This minimum occurs for the smallest possible values of i and j that allow for an obstruction.
# The first case is an observer in row j=3, with an obstruction in row i=2.

critical_j = 3
critical_i = 2

# Calculate the radii for this critical configuration.
r_1 = initial_radius
r_i_crit = initial_radius + (critical_i - 1) * row_depth
r_j_crit = initial_radius + (critical_j - 1) * row_depth

# Calculate the upper bound for K from this critical case.
k_bound = 2 * (r_i_crit - r_1) * (r_j_crit - r_1)

# Step 5: The maximum integer value for K is the largest integer strictly less than k_bound.
# We can find this using math.ceil(k_bound) - 1. For example, if k_bound=9.0, max int K is 8. If k_bound=9.1, max int K is 9.
max_K = math.ceil(k_bound) - 1

# Print the final result, showing the numbers used in the equation as requested.
print("The problem requires finding the maximum integer K for the parliament's parabolic design (h = r^2 / K).")
print("Visibility for all members is required, leading to the inequality: K < 2 * (r_i - r_1) * (r_j - r_1)")
print("where r_j is the observer's row radius, and r_i is the obstruction's row radius.\n")
print(f"There are {total_members} members in {num_sections} sections, which means there are {num_rows} rows.")
print("To ensure the inequality holds for all members, we must find the most restrictive (minimum) value for the right side of the equation.")
print(f"This critical case occurs for the first possible obstruction, with the observer at row j={critical_j} and the obstruction at row i={critical_i}.\n")
print("The values for the radii in this case are:")
print(f"  - Speaker's row radius (r_1): {r_1:.1f} m")
print(f"  - Obstruction row radius (r_2): {r_i_crit:.1f} m")
print(f"  - Observer's row radius (r_3): {r_j_crit:.1f} m\n")
print("Plugging these values into the inequality gives the final constraint:")
print(f"K < 2 * (r_2 - r_1) * (r_3 - r_1)")
print(f"K < 2 * ({r_i_crit:.1f} - {r_1:.1f}) * ({r_j_crit:.1f} - {r_1:.1f})")
print(f"K < 2 * ({r_i_crit - r_1:.1f}) * ({r_j_crit - r_1:.1f})")
print(f"K < {k_bound}\n")
print(f"Since K must be an integer, the maximum value it can take is {max_K}.")
print(f"<<<{max_K}>>>")
