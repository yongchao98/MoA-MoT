import math

# Billet dimensions
billet_l = 8
billet_w = 8
billet_h = 4

# Product specifications
# B2: Ball, 2cm radius
b2_radius = 2
b2_price = 150
b2_block_side = b2_radius * 2

# B1: Ball, 1cm diameter
b1_diameter = 1
b1_price = 1

# T1: Cube, 0.8cm side
t1_side = 0.8
t1_price = 1

# --- Step 1 & 2: Calculate maximum B2 production ---
# The most valuable product is B2. We prioritize it.
# Calculate how many 4x4x4 blocks can fit in the 8x8x4 billet.
num_b2 = math.floor(billet_l / b2_block_side) * math.floor(billet_w / b2_block_side) * math.floor(billet_h / b2_block_side)
value_from_b2 = num_b2 * b2_price

# --- Step 3: Calculate waste material ---
# Volume of the blocks used for B2s
total_block_volume_for_b2 = num_b2 * (b2_block_side ** 3)
# Volume of the actual B2 balls
volume_of_one_b2 = (4/3) * math.pi * (b2_radius ** 3)
total_ball_volume = num_b2 * volume_of_one_b2
# Total waste volume
total_waste_volume = total_block_volume_for_b2 - total_ball_volume

# --- Step 4 & 5: Maximize value from waste ---
# Determine which smaller product is more efficient in terms of value per volume.
# Value per cm^3 for B1: price / volume
b1_radius = b1_diameter / 2
b1_volume = (4/3) * math.pi * (b1_radius ** 3)
b1_value_density = b1_price / b1_volume # ~1.91

# Value per cm^3 for T1: price / volume
t1_volume = t1_side ** 3
t1_value_density = t1_price / t1_volume # ~1.95
# T1 is more valuable per unit of volume, so we use the waste to make T1 cubes.

# Calculate how many T1 cubes can be made from the waste.
num_t1_from_waste = math.floor(total_waste_volume / t1_volume)
value_from_t1 = num_t1_from_waste * t1_price

# --- Step 6: Calculate final total value ---
total_max_value = value_from_b2 + value_from_t1

# --- Print the final equation and result ---
# The equation shows the value from 4 B2 balls plus the value from T1 cubes made from the waste.
print(f"The optimal strategy is to cut {num_b2} B2 balls and use the waste to make T1 cubes.")
print("The final calculation for the maximum value is:")
print(f"{num_b2} * {b2_price} + floor(({num_b2} * {b2_block_side}^3 - {num_b2} * (4/3)*pi*{b2_radius}^3) / {t1_side}^3) * {t1_price} = {value_from_b2} + {num_t1_from_waste} = {total_max_value}")
