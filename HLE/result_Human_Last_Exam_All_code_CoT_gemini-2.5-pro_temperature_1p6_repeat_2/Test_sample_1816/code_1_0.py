import math

# Step 1: Define initial parameters and calculate the volume of a full cube.
side_length = 10  # in mm
volume_of_one_cube = side_length ** 3

print(f"The side length of a solid aluminum cube is {side_length} mm.")
print(f"The volume of one full cube is {side_length} * {side_length} * {side_length} = {volume_of_one_cube} mm^3.\n")

# Step 2: Calculate the volume of material removed from a single chamfered edge.
# The chamfer has a sloping width of sqrt(2) at 45 degrees.
# This means the cross-section is a right-angled isosceles triangle.
# Let 'a' be the length of the two equal sides of the triangle.
# By Pythagorean theorem: a^2 + a^2 = (sqrt(2))^2 => 2a^2 = 2 => a = 1 mm.
chamfer_cut_side = 1  # in mm
triangle_area = 0.5 * chamfer_cut_side * chamfer_cut_side
volume_per_chamfer = triangle_area * side_length

print(f"The material removed from one edge is a triangular prism.")
print(f"Its cross-section is a right triangle with two sides of length {chamfer_cut_side} mm.")
print(f"The volume of material removed per edge is (1/2 * {chamfer_cut_side} * {chamfer_cut_side}) * {side_length} = {volume_per_chamfer} mm^3.\n")

# Step 3: Calculate the total recycled material from one chamfered cube.
num_chamfered_edges = 4
total_recycled_volume_per_cube = num_chamfered_edges * volume_per_chamfer

print(f"One cube has {num_chamfered_edges} edges chamfered.")
print(f"The total recycled volume from one cube is {num_chamfered_edges} * {volume_per_chamfer} = {total_recycled_volume_per_cube} mm^3.\n")

# Step 4: Calculate the number of chamfered cubes needed.
# We need to find 'N' such that N * total_recycled_volume_per_cube = volume_of_one_cube
num_cubes_needed = volume_of_one_cube / total_recycled_volume_per_cube

print("To find how many chamfered cubes (N) are needed for one new cube, we solve the equation:")
print(f"N * {total_recycled_volume_per_cube} = {volume_of_one_cube}")
print(f"N = {volume_of_one_cube} / {total_recycled_volume_per_cube}")
print(f"N = {int(num_cubes_needed)}\n")
print(f"Therefore, {int(num_cubes_needed)} chamfered cubes are needed.")

<<<50>>>