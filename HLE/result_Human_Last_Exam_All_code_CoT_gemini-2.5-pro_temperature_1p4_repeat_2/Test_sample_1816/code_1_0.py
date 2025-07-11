import math

# Step 1: Define initial parameters and calculate the volume of the original cube.
side_length = 10  # in mm
volume_cube = side_length ** 3
print(f"The side length of the cube is {side_length} mm.")
print(f"The volume of one original cube is {side_length}^3 = {volume_cube} mm^3.\n")

# Step 2: Calculate the dimensions of the cut-off material for one chamfer.
# The chamfer has a sloping edge (hypotenuse of the cross-section triangle) of sqrt(2).
# The cut is at 45 degrees, making the cross-section an isosceles right-angled triangle.
# Let 'x' be the length of the two equal sides of the triangle.
# By Pythagorean theorem: x^2 + x^2 = (sqrt(2))^2
# 2 * x^2 = 2  => x^2 = 1 => x = 1 mm.
chamfer_leg_length = 1  # in mm
print("A chamfer removes a triangular prism of material from an edge.")
print(f"The cross-section of this prism is a right-angled triangle with sides of length {chamfer_leg_length} mm, {chamfer_leg_length} mm, and sqrt(2) mm.\n")

# Step 3: Calculate the volume of one chamfer cut.
# The area of the triangular cross-section is (1/2) * base * height.
cross_section_area = 0.5 * chamfer_leg_length * chamfer_leg_length
# The volume of the removed prism is the cross-section area times the length of the edge.
volume_one_chamfer = cross_section_area * side_length
print(f"The volume of material removed from one chamfer is:")
print(f"Volume = (0.5 * {chamfer_leg_length} * {chamfer_leg_length}) * {side_length} = {volume_one_chamfer} mm^3.\n")

# Step 4: Calculate the total recycled volume from one cube.
num_chamfered_edges = 4
total_recycled_volume_per_cube = num_chamfered_edges * volume_one_chamfer
print(f"One face has {num_chamfered_edges} edges chamfered.")
print(f"Total recycled volume per cube = {num_chamfered_edges} * {volume_one_chamfer} = {total_recycled_volume_per_cube} mm^3.\n")

# Step 5: Calculate the number of cubes needed.
# This is the volume of a full cube divided by the recycled volume per cube.
# We use math.ceil to round up to the nearest whole number, as we need enough material.
num_cubes_needed = math.ceil(volume_cube / total_recycled_volume_per_cube)

print("To find the number of chamfered cubes needed to make a new cube, we use the following equation:")
print(f"Number of Cubes = (Volume of one cube) / (Total recycled volume per cube)")
print(f"Number of Cubes = {volume_cube} / ({num_chamfered_edges} * {volume_one_chamfer})")
print(f"Number of Cubes = {volume_cube} / {total_recycled_volume_per_cube}")
print(f"Calculated result: {num_cubes_needed}\n")

print(f"Therefore, {int(num_cubes_needed)} chamfered cubes are needed to manufacture another cube.")
<<<50>>>