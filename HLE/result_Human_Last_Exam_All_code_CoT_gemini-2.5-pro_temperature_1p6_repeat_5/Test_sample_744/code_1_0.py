import math

# 1. Define constants for the material and the product
tube_outer_dim = 20  # cm
tube_thickness = 4   # cm
tube_length = 100    # cm (converted from 1m)
ball_radius = 2      # cm

# Values for the products
value_whole_ball = 3
value_welded_ball = 2

# 2. Derive dimensions based on the given values
tube_inner_dim = tube_outer_dim - 2 * tube_thickness
ball_diameter = 2 * ball_radius

# We assume each ball is machined from a cube of material with a side equal to the ball's diameter.
block_side = ball_diameter

# 3. Calculate how many blocks can be cut from the tube

# First, find the number of blocks in the tube's cross-section.
# We can model this as a grid.
outer_grid_side_count = tube_outer_dim // block_side
inner_grid_side_count = tube_inner_dim // block_side

# The number of blocks in the cross-section is the outer area minus the inner (hollow) area.
num_blocks_in_cross_section = (outer_grid_side_count ** 2) - (inner_grid_side_count ** 2)

# Next, find the number of layers of blocks along the tube's length.
num_layers = tube_length // block_side

# Finally, calculate the total number of blocks (which is the total number of balls we can make).
total_balls = num_blocks_in_cross_section * num_layers

# 4. Determine the optimal strategy and calculate the maximum value.
# Since a whole ball (value=3) is more valuable than a welded ball (value=2)
# for the same amount of material, the optimal strategy is to only make whole balls.
max_value = total_balls * value_whole_ball

# 5. Print the step-by-step thinking process and the final calculation.
print("--- Step 1: Calculate the number of balls that fit in the cross-section ---")
print(f"The number of {int(block_side)}x{int(block_side)}cm material blocks in the cross-section is calculated as:")
print(f"({tube_outer_dim}/{int(block_side)})^2 - ({tube_inner_dim}/{int(block_side)})^2 = {outer_grid_side_count}^2 - {inner_grid_side_count}^2 = {num_blocks_in_cross_section}")
print("")
print("--- Step 2: Calculate the number of layers along the length ---")
print(f"The number of {int(block_side)}cm layers along the {tube_length}cm length is:")
print(f"{tube_length} / {int(block_side)} = {num_layers}")
print("")
print("--- Step 3: Calculate the total number of manufacturable balls ---")
print(f"Total balls = (Blocks in cross-section) * (Number of layers)")
print(f"Total balls = {num_blocks_in_cross_section} * {num_layers} = {total_balls}")
print("")
print("--- Step 4: Calculate the maximum value ---")
print("To maximize value, we produce only whole balls (value=3) instead of welded balls (value=2).")
print("The final calculation is:")
print(f"{total_balls} * {value_whole_ball} = {max_value}")
