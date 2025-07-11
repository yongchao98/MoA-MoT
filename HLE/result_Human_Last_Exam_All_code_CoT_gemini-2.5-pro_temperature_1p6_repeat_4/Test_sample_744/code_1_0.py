import math

# Step 1: Define the dimensions of the tube and the balls
outer_side = 20.0  # cm
length = 100.0  # cm (1m)
thickness = 4.0  # cm
ball_radius = 2.0  # cm

# Calculate derived dimensions
ball_diameter = ball_radius * 2
inner_side = outer_side - 2 * thickness

print(f"Hollow Tube Dimensions:")
print(f"  - Outer side: {outer_side} cm")
print(f"  - Inner side: {inner_side} cm")
print(f"  - Wall thickness: {thickness} cm")
print(f"  - Length: {length} cm")
print("\nBall Dimensions:")
print(f"  - Radius: {ball_radius} cm")
print(f"  - Diameter: {ball_diameter} cm")

# Step 2: Determine how many whole balls can be packed into the material
# The ball diameter is 4cm, which perfectly matches the wall thickness.
# We can think of the cross-section of the tube as a grid of 4x4cm cells.

# Number of 4x4 cells in a 20x20 cross-section
cells_in_outer_square = (outer_side / ball_diameter)**2

# Number of 4x4 cells in the 12x12 inner hole
cells_in_inner_square = (inner_side / ball_diameter)**2

# Number of 4x4 cells in the material's cross-section
cells_per_slice = cells_in_outer_square - cells_in_inner_square

# Number of 4cm thick slices along the 100cm length
num_slices = length / ball_diameter

# Total number of 4x4x4cm cubes that can be cut, each yielding one whole ball
total_whole_balls = cells_per_slice * num_slices

print("\nPacking Calculation:")
print(f"A {outer_side}x{outer_side} cm cross-section can fit {cells_in_outer_square:.0f} balls across.")
print(f"The {inner_side}x{inner_side} cm hole removes {cells_in_inner_square:.0f} of those ball locations.")
print(f"Number of balls per 4cm-thick slice: {cells_per_slice:.0f}")
print(f"Number of 4cm slices along the {length}cm length: {num_slices:.0f}")
print(f"Total number of whole balls that can be made: {total_whole_balls:.0f}")

# Step 3: Evaluate the production options and calculate the maximum value
value_whole_ball = 3
value_welded_ball = 2

# To make one whole ball, we need a 4x4x4 cm block of material. This yields a value of 3.
# To make one welded ball (from two half-balls), we also need a 4x4x4 cm block of material. This yields a value of 2.
# Since 3 > 2, it is always better to make a whole ball.

# Step 4: Calculate the maximum value
max_value = total_whole_balls * value_whole_ball

print("\nValue Calculation:")
print("Since making a whole ball (value 3) is more valuable than making a welded ball (value 2) from the same amount of material, we should only make whole balls.")
print(f"The maximum value is calculated by multiplying the total number of whole balls by their value.")
print("\nFinal Equation:")
print(f"{int(total_whole_balls)} * {value_whole_ball} = {int(max_value)}")

<<<C>>>