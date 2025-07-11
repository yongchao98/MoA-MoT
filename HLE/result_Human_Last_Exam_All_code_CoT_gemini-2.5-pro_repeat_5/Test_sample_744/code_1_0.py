import math

# --- Problem Parameters ---
# Tube dimensions
outer_dim = 20  # cm
length = 100    # cm (1m)
thickness = 4   # cm

# Ball dimensions
ball_radius = 2 # cm
ball_diameter = ball_radius * 2

# Value of each type of ball
whole_ball_value = 3
welded_ball_value = 2

# --- Step 1: Calculate material geometry ---
# We determine the volume of material by considering how many
# bounding boxes of the balls can be cut from it.
print("Step 1: Calculate material and product dimensions.")
inner_dim = outer_dim - 2 * thickness
material_cross_section_area = outer_dim**2 - inner_dim**2
print(f"The material's cross-sectional area is ({outer_dim}x{outer_dim}) - ({inner_dim}x{inner_dim}) = {material_cross_section_area} cm^2.")

# --- Step 2: Calculate how many balls can be cut ---
# The space needed to cut a ball is a cube with side length equal to the ball's diameter.
ball_bounding_box_area = ball_diameter**2
print(f"A ball has a diameter of {ball_diameter} cm and requires a 4x4 cm area in the cross-section.")
print("-" * 30)

print("Step 2: Calculate the number of balls that can be manufactured.")
# Calculate how many rows of balls can be placed in the cross-section.
# This is the material's cross-section area divided by the ball's bounding box area.
num_rows = material_cross_section_area // ball_bounding_box_area
print(f"Number of ball rows fitting in the cross-section: {material_cross_section_area} / {ball_bounding_box_area} = {num_rows}")

# Calculate how many balls can fit along the length of the tube.
num_balls_per_row = length // ball_diameter
print(f"Number of balls fitting along the 100cm length: {length} / {ball_diameter} = {num_balls_per_row}")
print("-" * 30)

# --- Step 3: Determine the optimal manufacturing strategy ---
print("Step 3: Determine the optimal strategy.")
# The total number of whole balls we can make.
total_possible_whole_balls = num_rows * num_balls_per_row
print(f"Total whole balls possible: {num_rows} * {num_balls_per_row} = {total_possible_whole_balls}")

# Because the dimensions of the material are perfectly divisible by the ball's diameter
# (both in cross-section and length), there is no leftover material.
# Therefore, we cannot make any extra half-balls.
# The choice is to make a whole ball (value 3) or a welded ball (value 2) from each available slot.
# To maximize value, we should only make whole balls.
num_whole_balls_to_make = total_possible_whole_balls
num_welded_balls_to_make = 0
print("Since whole balls (value 3) are more valuable than welded balls (value 2), we will only make whole balls.")
print("-" * 30)

# --- Step 4: Calculate the maximum value ---
print("Step 4: Calculate the final maximum value.")
max_value = (num_whole_balls_to_make * whole_ball_value) + (num_welded_balls_to_make * welded_ball_value)

print("Final Equation:")
print(f"{num_whole_balls_to_make} * {whole_ball_value} = {max_value}")
