import math

# Step 1: Define the initial parameters
side_length = 10  # mm
chamfer_hypotenuse = math.sqrt(2)  # mm
num_chamfered_edges = 4

# Step 2: Calculate the dimensions of the triangular cross-section of the removed material
# For a 45-degree chamfer, the cross-section is an isosceles right-angled triangle.
# Let 'a' be the length of the two equal sides (the legs of the triangle).
# a^2 + a^2 = chamfer_hypotenuse^2
# 2 * a^2 = (sqrt(2))^2 = 2
# a^2 = 1 => a = 1 mm
leg_length_a = 1

# Step 3: Calculate the volume of one removed piece (a triangular prism)
# Area of the triangular base = 0.5 * base * height = 0.5 * a * a
area_triangle_base = 0.5 * leg_length_a * leg_length_a
# Volume = area * length (length is the side_length of the cube)
volume_one_chamfer = area_triangle_base * side_length

# Step 4: Calculate the total recycled volume from one chamfered cube
recycled_volume_per_cube = num_chamfered_edges * volume_one_chamfer

# Step 5: Calculate the volume of a new cube
volume_new_cube = side_length ** 3

# Step 6: Calculate how many chamfered cubes are needed
# This is the volume of a new cube divided by the recycled volume from one cube.
num_cubes_needed = volume_new_cube / recycled_volume_per_cube

# Print the final equation and the result
print("The volume of a new cube is {}^3 = {} mm^3.".format(side_length, int(volume_new_cube)))
print("The volume of recycled material from one chamfered cube is {} * (0.5 * {} * {} * {}) = {} mm^3.".format(num_chamfered_edges, leg_length_a, leg_length_a, side_length, int(recycled_volume_per_cube)))
print("\nTo find the number of chamfered cubes needed, we solve the equation:")
print("{} / {} = {}".format(int(volume_new_cube), int(recycled_volume_per_cube), int(num_cubes_needed)))
print("\nTherefore, {} chamfered cubes are needed to manufacture another cube from the recycled material.".format(int(num_cubes_needed)))
