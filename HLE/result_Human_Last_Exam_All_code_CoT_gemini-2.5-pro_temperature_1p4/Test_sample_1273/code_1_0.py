import math

# Step 1: Define billet and product specifications
billet_l, billet_w, billet_h = 8, 8, 4
billet_volume = billet_l * billet_w * billet_h

# Product B2: Ball, radius 2cm, price 150
b2_price = 150
b2_radius = 2
b2_diameter = b2_radius * 2
b2_volume = (4/3) * math.pi * (b2_radius**3)

# Product T1: Cube, side 0.8cm, price 1
t1_price = 1
t1_side = 0.8
t1_volume = t1_side**3

# Product B1: Ball, diameter 1cm, price 1
b1_price = 1
b1_diameter = 1
b1_radius = b1_diameter / 2
b1_volume = (4/3) * math.pi * (b1_radius**3)

# Step 2: Prioritize cutting the most valuable item, B2 balls
# Calculate the maximum number of B2 balls that can be cut from the billet
num_b2_l = math.floor(billet_l / b2_diameter)
num_b2_w = math.floor(billet_w / b2_diameter)
num_b2_h = math.floor(billet_h / b2_diameter)
max_b2_cut = num_b2_l * num_b2_w * num_b2_h
value_from_b2 = max_b2_cut * b2_price

# Step 3: Calculate the volume of leftover steel after cutting B2 balls
volume_used_by_b2_balls = max_b2_cut * b2_volume
volume_of_waste_material = billet_volume - volume_used_by_b2_balls

# Step 4: Determine the best use of the waste material
# Compare value per unit volume for T1 and B1
# Price/Vol for T1 = 1 / 0.8^3 = 1.953
# Price/Vol for B1 = 1 / (4/3*pi*0.5^3) = 1.910
# T1 cubes are the more efficient use of the leftover steel.
num_t1_from_waste = math.floor(volume_of_waste_material / t1_volume)
value_from_t1 = num_t1_from_waste * t1_price

# Step 5: Calculate the final total value
total_value = value_from_b2 + value_from_t1

# Print the final equation with all the numbers
print("Total Value = (Number of B2 balls * Price of B2) + (Number of T1 cubes from waste * Price of T1)")
print(f"Total Value = ({max_b2_cut} * {b2_price}) + ({num_t1_from_waste} * {t1_price}) = {value_from_b2} + {value_from_t1} = {total_value}")
