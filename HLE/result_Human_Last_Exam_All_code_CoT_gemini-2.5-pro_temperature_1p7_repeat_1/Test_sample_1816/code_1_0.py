import math

# Step 1: Define initial parameters and calculate the volume of a standard cube.
cube_side_length = 10  # in mm
volume_of_one_cube = cube_side_length ** 3
print(f"The volume of one standard cube is {cube_side_length}^3 = {volume_of_one_cube} mm^3.")

# Step 2: Calculate the volume of material removed by one chamfer.
# The cross-section of the cut is a right-angled triangle.
# The hypotenuse is sqrt(2). Let the other two equal sides be 'a'.
# a^2 + a^2 = (sqrt(2))^2  => 2*a^2 = 2 => a = 1 mm.
cut_triangle_side = 1 # in mm
cut_triangle_area = 0.5 * cut_triangle_side * cut_triangle_side
# The length of the cut piece is the side length of the cube.
volume_of_one_chamfer = cut_triangle_area * cube_side_length
print(f"The volume of material removed from one chamfered edge is 0.5 * {cut_triangle_side} * {cut_triangle_side} * {cube_side_length} = {volume_of_one_chamfer} mm^3.")

# Step 3: Calculate the total recycled volume from one chamfered cube.
num_chamfered_edges = 4
recycled_volume_per_cube = num_chamfered_edges * volume_of_one_chamfer
print(f"A total of {num_chamfered_edges} edges are chamfered on each cube.")
print(f"The total recycled volume per cube is {num_chamfered_edges} * {volume_of_one_chamfer} = {recycled_volume_per_cube} mm^3.")

# Step 4: Calculate how many chamfered cubes are needed.
# This is the volume of a new cube divided by the recycled volume from one cube.
num_cubes_needed = volume_of_one_cube / recycled_volume_per_cube
print("\nTo find the number of chamfered cubes needed to make one new cube, we perform the following calculation:")
print(f"Number of cubes = (Volume of one cube) / (Recycled volume per cube)")
print(f"Number of cubes = {volume_of_one_cube} / {recycled_volume_per_cube} = {int(num_cubes_needed)}")
