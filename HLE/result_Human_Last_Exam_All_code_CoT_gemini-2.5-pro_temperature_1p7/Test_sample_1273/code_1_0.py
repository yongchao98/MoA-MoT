import math

# Billet dimensions
billet_l, billet_w, billet_h = 8, 8, 4

# Product B2: Ball, radius 2cm, price 150
b2_radius = 2
b2_diameter = b2_radius * 2
b2_price = 150

# Product T1: Cube, side 0.8cm, price 1
t1_side = 0.8
t1_price = 1

# Step 1: Calculate how many B2 balls can be cut from the billet
# We need to fit the bounding box (a cube of side = diameter) for each ball.
num_b2_l = math.floor(billet_l / b2_diameter)
num_b2_w = math.floor(billet_w / b2_diameter)
num_b2_h = math.floor(billet_h / b2_diameter)
num_b2 = num_b2_l * num_b2_w * num_b2_h

# Calculate value from B2 balls
value_from_b2 = num_b2 * b2_price

# Step 2: Calculate the remaining material volume
billet_volume = billet_l * billet_w * billet_h
b2_sphere_volume = (4/3) * math.pi * (b2_radius ** 3)
total_b2_material_volume = num_b2 * b2_sphere_volume
remaining_volume = billet_volume - total_b2_material_volume

# Step 3: Use the remaining volume to make the most efficient smaller product (T1)
t1_cube_volume = t1_side ** 3
num_t1 = math.floor(remaining_volume / t1_cube_volume)

# Calculate value from T1 cubes
value_from_t1 = num_t1 * t1_price

# Step 4: Calculate the total combined value
total_value = value_from_b2 + value_from_t1

# Print the final equation as requested
print(f"The optimal combination found is {num_b2} B2 balls and {num_t1} T1 cubes.")
print("The equation for the total value is:")
print(f"{num_b2} * {b2_price} + {num_t1} * {t1_price} = {total_value}")