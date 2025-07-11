import math

# --- Define Billet and Product Properties ---

# Billet dimensions (cm)
billet_l, billet_w, billet_h = 8, 8, 4

# B2 Ball properties
b2_price = 150
b2_radius = 2
b2_block_side = b2_radius * 2  # B2 ball is cut from a 4x4x4 cm cube

# T1 Cube properties
t1_price = 1
t1_side = 0.8

# --- Step 1: Calculate the number of primary products (B2 balls) ---

# We can fit floor(8/4)xfloor(8/4)xfloor(4/4) = 2x2x1 = 4 B2 balls.
num_b2 = (billet_l // b2_block_side) * (billet_w // b2_block_side) * (billet_h // b2_block_side)

# --- Step 2: Calculate the waste volume after cutting B2 balls ---

# The 4 blocks for the B2 balls occupy the entire billet.
billet_volume = billet_l * billet_w * billet_h

# Volume of the material that becomes B2 balls
b2_ball_volume = (4/3) * math.pi * (b2_radius ** 3)
total_b2_volume = num_b2 * b2_ball_volume

# The waste is the volume of the billet minus the volume of the spheres
waste_volume = billet_volume - total_b2_volume

# --- Step 3: Calculate how many T1 cubes can be made from waste ---

# Volume of a single T1 cube
t1_cube_volume = t1_side ** 3

# Number of T1 cubes from waste (assuming waste can be fully utilized)
num_t1 = math.floor(waste_volume / t1_cube_volume)

# --- Step 4: Calculate the total maximized value ---

# Total value = (Value from B2 balls) + (Value from T1 cubes)
total_value = (num_b2 * b2_price) + (num_t1 * t1_price)

# --- Print the final equation and result ---

print("To maximize the value, we first cut the most profitable items, B2 balls, and then use the leftover material to make T1 cubes.")
print("The calculation for the total value is as follows:")
print(f"({num_b2} * {b2_price}) + ({num_t1} * {t1_price}) = {int(total_value)}")