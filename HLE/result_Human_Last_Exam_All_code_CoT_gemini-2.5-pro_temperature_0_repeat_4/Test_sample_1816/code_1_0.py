import math

# Define the properties of the cube and the chamfer process
side_length = 10  # in mm
num_chamfered_edges = 4
chamfer_hypotenuse_width = math.sqrt(2)  # in mm

# Step 1: Calculate the volume of a single, original cube
volume_cube = side_length**3

# Step 2: Calculate the volume of material removed by one chamfer
# The cross-section of the chamfer is a right-angled isosceles triangle.
# The chamfer width is the hypotenuse. Let 'a' be the length of the two equal sides.
# a^2 + a^2 = chamfer_hypotenuse_width^2  => 2*a^2 = 2 => a^2 = 1 => a = 1
cut_depth = math.sqrt(chamfer_hypotenuse_width**2 / 2)

# The volume of one removed piece (a triangular prism) is the area of the
# triangular base times the length of the edge.
volume_one_chamfer = (0.5 * cut_depth * cut_depth) * side_length

# Step 3: Calculate the total volume of recycled material from one cube
total_volume_removed_per_cube = num_chamfered_edges * volume_one_chamfer

# Step 4: Calculate how many chamfered cubes are needed to make one new cube
# This is the volume of a full cube divided by the recycled volume per cube.
num_cubes_needed = volume_cube / total_volume_removed_per_cube

# Print the results and the final equation
print(f"The volume of one original cube is {int(volume_cube)} mm^3.")
print(f"The volume of material removed from one chamfered edge is {volume_one_chamfer:.1f} mm^3.")
print(f"The total volume of recycled material from one cube is {int(total_volume_removed_per_cube)} mm^3.")
print("\nTo find the number of cubes needed, we use the following equation:")
print("Number of Cubes = (Volume of one cube) / (Total volume removed per cube)")
print("\nThe final equation with the calculated values is:")
print(f"{int(num_cubes_needed)} = {int(volume_cube)} / ({int(num_chamfered_edges)} * {int(volume_one_chamfer)})")
print(f"\nTherefore, {int(num_cubes_needed)} chamfered cubes are needed to recycle enough material for one new cube.")