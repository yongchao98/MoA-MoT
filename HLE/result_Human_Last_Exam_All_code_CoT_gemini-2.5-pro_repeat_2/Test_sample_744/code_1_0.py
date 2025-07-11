import math

# Step 1: Define the dimensions and values from the problem
# Tube dimensions in cm
tube_outer_side = 20
tube_length = 100  # 1m = 100cm
tube_thickness = 4

# Ball dimensions in cm
ball_radius = 2
ball_diameter = ball_radius * 2

# Values for the products
value_whole_ball = 3
value_welded_ball = 2

# Step 2: Calculate the dimensions of the inner hollow space
tube_inner_side = tube_outer_side - 2 * tube_thickness

# Step 3: Calculate how many balls can be cut from the material.
# The ball diameter (4cm) is the same as the tube's wall thickness.
# We can model the tube as four rectangular walls.

# Calculate how many balls fit along the length of the tube
balls_along_length = tube_length // ball_diameter

# Calculate balls from the two wider walls (e.g., top and bottom)
# These walls have a surface of 20cm x 100cm
balls_along_wide_wall_width = tube_outer_side // ball_diameter
balls_from_two_wide_walls = 2 * balls_along_wide_wall_width * balls_along_length

# Calculate balls from the two narrower walls (e.g., sides)
# These walls have a surface of 12cm x 100cm
balls_along_narrow_wall_width = tube_inner_side // ball_diameter
balls_from_two_narrow_walls = 2 * balls_along_narrow_wall_width * balls_along_length

# Calculate the total number of balls
total_balls = balls_from_two_wide_walls + balls_from_two_narrow_walls

# Step 4: Determine the best strategy and calculate the maximum value.
# For each 4x4x4cm block of material, we can either make:
# - A whole ball (value 3)
# - A welded ball (value 2)
# To maximize value, we should always choose to make whole balls.

max_value = total_balls * value_whole_ball

# Step 5: Print the final result and the equation used.
print("To maximize the value, we should only produce whole balls, as their value (3) is higher than welded balls (2).")
print("\nThe calculation for the total number of balls is:")
print(f"({balls_along_wide_wall_width} balls wide * {balls_along_length} balls long * 2 walls) + ({balls_along_narrow_wall_width} balls wide * {balls_along_length} balls long * 2 walls) = {total_balls} balls")

print("\nThe final calculation for the maximum value is:")
print(f"{total_balls} * {value_whole_ball} = {max_value}")
<<<C>>>