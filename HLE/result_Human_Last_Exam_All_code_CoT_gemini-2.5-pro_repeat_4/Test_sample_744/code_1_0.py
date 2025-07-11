import math

# Step 1: Define material and product dimensions and values
outer_side = 20  # cm
length = 100  # cm (1m)
thickness = 4  # cm
ball_radius = 2  # cm
value_whole_ball = 3
value_welded_ball = 2

print("--- Step 1: Analyzing the Dimensions ---")

# Step 2: Calculate derived dimensions
inner_side = outer_side - 2 * thickness
ball_diameter = 2 * ball_radius

print(f"The hollow tube has an outer side of {outer_side}cm, an inner side of {inner_side}cm, and a length of {length}cm.")
print(f"Each ball has a radius of {ball_radius}cm, which means a diameter of {ball_diameter}cm.")
print("To cut a whole ball, we need a 4x4x4 cm cube of material.")
print("-" * 20)

print("\n--- Step 2: Calculating the Number of Whole Balls ---")
# Step 3: Calculate how many balls can be cut from the material

# Calculate how many 4x4 squares fit in the cross-section
# We can model this as the number of cells in an outer grid minus the number of cells in the inner (hollow) grid.
outer_grid_dim = outer_side // ball_diameter
inner_grid_dim = inner_side // ball_diameter

balls_in_cross_section = (outer_grid_dim ** 2) - (inner_grid_dim ** 2)

# Calculate how many 4cm layers fit along the length
layers_along_length = length // ball_diameter

# Calculate total number of whole balls
total_whole_balls = balls_in_cross_section * layers_along_length

print(f"The number of 4cm sections that fit across the outer side is: {outer_side} / {ball_diameter} = {outer_grid_dim}")
print(f"The number of 4cm sections that fit across the inner side is: {inner_side} / {ball_diameter} = {inner_grid_dim}")
print(f"So, the number of balls that can be placed in a single cross-sectional slice is: ({outer_grid_dim}^2) - ({inner_grid_dim}^2) = {balls_in_cross_section}")
print(f"The number of 4cm layers that fit along the 100cm length is: {length} / {ball_diameter} = {layers_along_length}")
print(f"Total number of whole balls = {balls_in_cross_section} * {layers_along_length} = {total_whole_balls}")
print("-" * 20)


print("\n--- Step 3: Checking for Leftover Material for Welded Balls ---")
# Step 4: Check for leftover material
# Since the dimensions of the tube (100, 20, 12) are all perfectly divisible by the ball diameter (4),
# the entire volume is used up by the 4x4x4 cutting cubes for the whole balls.
# There is no leftover material to make half-balls.
total_welded_balls = 0
print("The material dimensions are perfectly divisible by the ball diameter.")
print(f"Therefore, there is no leftover material to create half-balls. Number of welded balls = {total_welded_balls}")
print("-" * 20)


print("\n--- Step 4: Calculating the Final Value ---")
# Step 5: Calculate the maximum value
max_value = (total_whole_balls * value_whole_ball) + (total_welded_balls * value_welded_ball)

print(f"The highest value is calculated by maximizing the number of more valuable whole balls.")
print(f"Total Value = (Number of Whole Balls * Value of Whole Ball) + (Number of Welded Balls * Value of Welded Ball)")
print(f"Total Value = ({total_whole_balls} * {value_whole_ball}) + ({total_welded_balls} * {value_welded_ball}) = {max_value}")

<<<C>>>