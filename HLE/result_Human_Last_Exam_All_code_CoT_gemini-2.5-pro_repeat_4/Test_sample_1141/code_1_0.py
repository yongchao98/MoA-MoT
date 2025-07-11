# Define prices for each product type
price_b2 = 150
price_b1 = 1
price_t1 = 5

# --- Scenario A: T1 cubes only ---
# The non-overlapping constraint for T1 cubes allows them to be placed on a grid.
# The centers must be separated by at least 2 units on each axis.
# We can place them at coordinates (1, 3, 5, ...), (1, 3, 5, ...), (1, 3, 5, ...).

# T1 center coordinate ranges
x_range_t1 = (1, 31)
y_range_t1 = (1, 21)
z_range_t1 = (1, 7)

# Step size for the grid
step = 2

# Calculate the number of possible positions on each axis
num_x_t1 = (x_range_t1[1] - x_range_t1[0]) // step + 1
num_y_t1 = (y_range_t1[1] - y_range_t1[0]) // step + 1
num_z_t1 = (z_range_t1[1] - z_range_t1[0]) // step + 1

# Total number of T1 cubes
num_t1 = num_x_t1 * num_y_t1 * num_z_t1

# As determined in the analysis, if we place this many T1s, no space is left for B1s
# according to the given constraints. Also, no B2s can be placed with T1s.
num_b2_scenario_A = 0
num_b1_scenario_A = 0

# Calculate total value for this scenario
total_value_A = (num_t1 * price_t1) + (num_b1_scenario_A * price_b1) + (num_b2_scenario_A * price_b2)

# --- Scenario B: B2 balls only ---
# From our analysis, this scenario is clearly suboptimal, but we include it for completeness.
# A simple 3x2 grid packing allows for 6 B2 balls.
num_b2_scenario_B = 6
# Even if we filled the rest of the billet with B1 balls, the value would be less than Scenario A.
# The maximum number of B1 balls is less than 704.
# Max possible value is < 6*150 + 704*1 = 900 + 704 = 1604, which is less than 3520.

# --- Conclusion ---
# Scenario A gives the highest value.
# The optimal solution is to cut only T1 cubes.

final_n_b2 = num_b2_scenario_A
final_n_b1 = num_b1_scenario_A
final_n_t1 = num_t1
max_value = total_value_A

# Print the final answer in the required format
print("Based on the provided formulation, the optimal strategy is to produce:")
print(f"- {final_n_b2} B2 balls")
print(f"- {final_n_b1} B1 balls")
print(f"- {final_n_t1} T1 cubes")
print("\nThe calculation for the maximum total value is:")
print(f"{final_n_b2} * {price_b2} + {final_n_b1} * {price_b1} + {final_n_t1} * {price_t1} = {max_value}")

<<<3520>>>