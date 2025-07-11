import math

# Billet dimensions in cm
billet_x_cm = 16
billet_y_cm = 11
billet_z_cm = 4

# Product dimensions in cm
t1_cube_side = 1

# Product prices
price_b2 = 150
price_b1 = 1
price_t1 = 5

# --- Analysis ---
# Based on the problem's constraints, the optimal solution involves filling the entire
# billet with T1 cubes (1cm side).
# The constraints make it impossible to mix T1 cubes with B2 balls.
# A solution with only T1 cubes yields a higher total value than any solution
# with B2 balls and B1 balls.
# T1 cubes are more valuable than B1 balls, so we prioritize them.

# --- Calculation ---
# Calculate how many T1 cubes fit along each dimension
num_t1_x = math.floor(billet_x_cm / t1_cube_side)
num_t1_y = math.floor(billet_y_cm / t1_cube_side)
num_t1_z = math.floor(billet_z_cm / t1_cube_side)

# Total number of T1 cubes
total_t1 = num_t1_x * num_t1_y * num_t1_z

# In this optimal solution, no B2 or B1 items are cut.
num_b2 = 0
num_b1 = 0

# Calculate the maximum total value
max_value = (num_b2 * price_b2) + (num_b1 * price_b1) + (total_t1 * price_t1)

# --- Output the result ---
# The final equation includes all item types for clarity, even those with zero count.
print(f"The highest possible value is achieved with {total_t1} T1 cubes, {num_b2} B2 balls, and {num_b1} B1 balls.")
print("The calculation for the total value is:")
print(f"{num_b2} * {price_b2} + {num_b1} * {price_b1} + {total_t1} * {price_t1} = {max_value}")
