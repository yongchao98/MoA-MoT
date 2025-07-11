import math

# Define the properties of the cubes and the shift
cube_side_length = 1.0
shift_per_axis = 0.5

# Define the boundaries of the two cubes
# Cube 1: from 0 to 1 on each axis
cube1_min = 0.0
cube1_max = 1.0

# Cube 2: shifted by 0.5 on each axis
cube2_min = shift_per_axis
cube2_max = cube_side_length + shift_per_axis

# Find the boundaries of the intersection volume
# The start of the intersection is the maximum of the two minimums
# The end of the intersection is the minimum of the two maximums
intersect_min = max(cube1_min, cube2_min)
intersect_max = min(cube1_max, cube2_max)

# Calculate the side length of the intersection cube
intersect_side_length = intersect_max - intersect_min

# Calculate the volume of the intersection
volume = intersect_side_length ** 3

# --- Output the results step-by-step ---

print("Step 1: Determine the geometry of the intersection.")
print(f"Let Cube 1 be at [0, 1] on each axis.")
print(f"A shift of {shift_per_axis} along the main diagonal places Cube 2 at [{shift_per_axis}, {shift_per_axis+cube_side_length}] on each axis.")
print("The intersection is a smaller cube defined by the overlap on each axis.")
print("")

print("Step 2: Calculate the side length of the intersection cube.")
print(f"The overlap on any axis is from {intersect_min} to {intersect_max}.")
side_x = intersect_side_length
side_y = intersect_side_length
side_z = intersect_side_length
print(f"The side length of the intersection cube is: {intersect_max} - {intersect_min} = {intersect_side_length}")
print("")

print("Step 3: Calculate the volume of the intersection.")
print("The volume is the product of its side lengths.")
print(f"Volume = {side_x} * {side_y} * {side_z}")
print(f"Volume = {volume}")
