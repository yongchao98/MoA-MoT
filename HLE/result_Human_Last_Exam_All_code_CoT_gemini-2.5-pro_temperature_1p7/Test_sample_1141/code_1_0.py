# Item information
price_b2 = 150
price_b1 = 1
price_t1 = 5

# --- Strategy 1: Calculate value from B2 balls only ---
# B2 (ball, 2cm radius) centers must be >= 8 units apart.
# Billet dimensions in 0.5cm units: 32x22x8
# B2 center ranges: x in [4, 28], y in [4, 18], z = 4
b2_x_range = 28 - 4
b2_y_range = 18 - 4
b2_separation = 8
num_b2_x = (b2_x_range // b2_separation) + 1
num_b2_y = (b2_y_range // b2_separation) + 1
total_b2 = num_b2_x * num_b2_y
value_from_b2 = total_b2 * price_b2

# --- Strategy 2: Calculate value from T1 cubes only ---
# T1 (cube, 1cm side) are 2x2x2 units.
# With the corrected standard packing constraint, their centers can be 2 units apart.
# T1 center ranges: x in [1, 31], y in [1, 21], z in [1, 7]
t1_x_range = 31 - 1
t1_y_range = 21 - 1
t1_z_range = 7 - 1
t1_separation = 2
num_t1_x = (t1_x_range // t1_separation) + 1
num_t1_y = (t1_y_range // t1_separation) + 1
num_t1_z = (t1_z_range // t1_separation) + 1
total_t1 = num_t1_x * num_t1_y * num_t1_z
value_from_t1 = total_t1 * price_t1

# --- Determine the optimal strategy ---
if value_from_t1 > value_from_b2:
    optimal_num_t1 = total_t1
    optimal_num_b2 = 0
    optimal_num_b1 = 0
    total_value = value_from_t1
else:
    # This part is for completeness, but T1 strategy is superior
    optimal_num_t1 = 0
    optimal_num_b2 = total_b2
    optimal_num_b1 = 0
    total_value = value_from_b2
    
# Print the final equation for the optimal solution
print(f"{optimal_num_t1} * {price_t1} + {optimal_num_b2} * {price_b2} + {optimal_num_b1} * {price_b1} = {total_value}")
