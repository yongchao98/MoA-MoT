import math

# --- Step 1: Define Dimensions ---
tube_outer_side_cm = 20
tube_length_cm = 100
tube_thickness_cm = 4
ball_radius_cm = 2

# --- Derived Dimensions ---
ball_diameter_cm = ball_radius_cm * 2
# The inner side is the outer side minus the thickness on both sides
tube_inner_side_cm = tube_outer_side_cm - 2 * tube_thickness_cm

# --- Step 2: Calculate the number of balls that can be made ---
# Since the material dimensions are perfect multiples of the ball diameter,
# we can determine the number of 4x4x4cm cubes we can cut.

# Calculate the number of ball 'slots' in the cross-section
material_cross_section_area = tube_outer_side_cm**2 - tube_inner_side_cm**2
ball_slot_cross_section_area = ball_diameter_cm**2
balls_per_layer = material_cross_section_area / ball_slot_cross_section_area

# Calculate the number of layers along the tube's length
num_layers = tube_length_cm / ball_diameter_cm

# Calculate the total number of whole balls
total_whole_balls = balls_per_layer * num_layers

# --- Step 3: Calculate the maximum value ---
# Value of a whole ball is 3. Value of a welded ball is 2.
# Since a whole ball has a higher value and the material can be perfectly
# utilized to create whole balls, we only make whole balls.
value_per_whole_ball = 3
max_value = total_whole_balls * value_per_whole_ball

# --- Step 4: Print the results ---
print("Step-by-step calculation:")
print(f"The material's cross-sectional area is ({tube_outer_side_cm}^2 - {tube_inner_side_cm}^2) = {int(material_cross_section_area)} cm^2.")
print(f"A ball requires a {int(ball_diameter_cm)}x{int(ball_diameter_cm)} cross-section, which is {int(ball_slot_cross_section_area)} cm^2.")
print(f"Number of balls per layer = {int(material_cross_section_area)} / {int(ball_slot_cross_section_area)} = {int(balls_per_layer)} balls.")
print(f"Number of layers along the {tube_length_cm}cm length = {tube_length_cm} / {int(ball_diameter_cm)} = {int(num_layers)} layers.")
print("-" * 20)
print("Final value calculation:")
print(f"Total whole balls that can be made = {int(balls_per_layer)} balls/layer * {int(num_layers)} layers = {int(total_whole_balls)} balls.")
print(f"Maximum Value = {int(total_whole_balls)} balls * {value_per_whole_ball} value/ball = {int(max_value)}")

<<<C>>>