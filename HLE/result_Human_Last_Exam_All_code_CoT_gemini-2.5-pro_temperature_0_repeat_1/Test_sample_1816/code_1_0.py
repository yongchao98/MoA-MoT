import math

# Step 1: Define the properties of the cube and the chamfer.
side_length = 10  # in mm
num_chamfered_edges = 4

# The volume of the original cube.
volume_cube = side_length ** 3

# Step 2: Calculate the volume of the material removed from one cube.
# The chamfer creates a sloping face with a width of sqrt(2) mm at a 45-degree angle.
# This means the cross-section of the removed material is a right-angled isosceles triangle.
# The width of the chamfer is the hypotenuse of this triangle.
# Let 'a' be the length of the two equal sides of the triangle (the cut depth on each face).
# Using the Pythagorean theorem: a^2 + a^2 = (sqrt(2))^2
# 2 * a^2 = 2, which means a = 1 mm.
cut_depth = 1

# The volume of material removed from one edge is a triangular prism.
# Its volume is the area of the triangular base times the length of the edge.
volume_per_chamfer = (0.5 * cut_depth * cut_depth) * side_length

# The total volume removed is for four edges.
total_removed_volume = volume_per_chamfer * num_chamfered_edges

# Step 3: Calculate how many chamfered cubes are needed to make one new cube.
# This is the volume of a full cube divided by the recycled volume from one chamfered cube.
num_cubes_needed = volume_cube / total_removed_volume

# Final step: Output the numbers in the final equation as requested.
# The equation is: Volume of a Cube / Recycled Volume per Cube = Number of Cubes Needed
print(f"{int(volume_cube)} / {int(total_removed_volume)} = {int(num_cubes_needed)}")