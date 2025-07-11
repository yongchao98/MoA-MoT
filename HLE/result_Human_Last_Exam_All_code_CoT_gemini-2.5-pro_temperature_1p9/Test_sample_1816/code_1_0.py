import math

# Step 1: Define variables and calculate the volume of a full cube.
side_length = 10  # in mm
volume_cube = side_length**3

# Step 2: Calculate the volume removed by a single chamfer.
# The chamfer process creates a cut-off piece that is a triangular prism.
# The cross-section is an isosceles right-angled triangle because the chamfer is at 45 degrees.
# The hypotenuse of this triangle is the sloping edge, which is sqrt(2) mm wide.
# Let 'x' be the length of the two equal sides of the triangle.
# From the Pythagorean theorem: x^2 + x^2 = (sqrt(2))^2 => 2*x^2 = 2 => x = 1 mm.
cut_depth = 1  # in mm
triangle_area = 0.5 * cut_depth * cut_depth
volume_one_chamfer = triangle_area * side_length

# Step 3: Calculate the total recycled volume from one chamfered cube.
num_chamfered_edges = 4
recycled_volume_per_cube = num_chamfered_edges * volume_one_chamfer

# Step 4: Calculate how many cubes are needed.
# This is the total volume needed (one cube) divided by the volume recovered from each cube.
num_cubes_needed = volume_cube / recycled_volume_per_cube

# Print the final equation with all the numbers.
print(f"The volume of a full cube is {int(volume_cube)} mm^3.")
print(f"The recycled material from one chamfered cube is {int(recycled_volume_per_cube)} mm^3.")
print("To find how many chamfered cubes are needed to manufacture one new cube, the final equation is:")
print(f"Number of Cubes = (Volume of a full cube) / (Recycled volume per cube)")
print(f"Number of Cubes = {int(volume_cube)} / {int(recycled_volume_per_cube)} = {int(num_cubes_needed)}")