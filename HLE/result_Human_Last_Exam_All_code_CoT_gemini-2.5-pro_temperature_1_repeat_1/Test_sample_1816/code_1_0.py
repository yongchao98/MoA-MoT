import math

# Step 1: Define initial parameters and calculate the volume of a single cube.
side_length = 10  # in mm
volume_cube = side_length**3
print(f"The volume of a single cube with side length {side_length} mm is {volume_cube} mm^3.\n")

# Step 2: Calculate the volume of material removed by one chamfer.
# The chamfer's cross-section is an isosceles right-angled triangle.
# The hypotenuse (sloping edge) is sqrt(2) mm.
# Let 'a' be the length of the other two equal sides.
# From the Pythagorean theorem: a^2 + a^2 = (sqrt(2))^2 => 2*a^2 = 2 => a = 1 mm.
cut_depth = 1.0 # in mm
area_triangle = 0.5 * cut_depth * cut_depth
volume_one_chamfer = area_triangle * side_length
print(f"The volume of material removed by a single chamfer is {volume_one_chamfer} mm^3.\n")

# Step 3: Calculate the total recycled volume from one cube.
# Four edges on one face are chamfered.
num_chamfered_edges = 4
volume_recycled_per_cube = num_chamfered_edges * volume_one_chamfer
print(f"The total recycled material from one cube (with {num_chamfered_edges} chamfers) is {volume_recycled_per_cube} mm^3.\n")

# Step 4: Calculate the number of cubes needed to form a new cube.
# This is the volume of a full cube divided by the recycled volume per cube.
num_cubes_needed = volume_cube / volume_recycled_per_cube
print("To find the number of chamfered cubes needed, we use the following equation:")
print(f"Number of Cubes = (Volume of a Cube) / (Recycled Volume per Cube)")
print(f"Number of Cubes = {int(volume_cube)} / {int(volume_recycled_per_cube)} = {int(num_cubes_needed)}")

<<<50>>>