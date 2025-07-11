import math

# Define the properties of the cube and the chamfer process
side_length = 10  # in mm
num_chamfered_edges = 4
chamfer_hypotenuse = math.sqrt(2)  # width of the sloping edge in mm

# Calculate the total volume of one original cube
volume_of_one_cube = side_length ** 3

# Calculate the volume of material removed from one cube
# The cross-section of the chamfer is a right-angled isosceles triangle.
# We find the length of the triangle's other two sides 'a' using Pythagoras: a^2 + a^2 = c^2
# 2 * a^2 = (sqrt(2))^2  => 2 * a^2 = 2  => a = 1 mm
cut_depth = math.sqrt(chamfer_hypotenuse**2 / 2)

# The volume of material removed from a single edge (a triangular prism)
# Volume = Area of triangular base * length of the prism
volume_per_chamfer = (0.5 * cut_depth * cut_depth) * side_length

# Total volume recycled from one cube
total_recycled_volume_per_cube = volume_per_chamfer * num_chamfered_edges

# Calculate the number of cubes needed to make one new cube from recycled material
# This is the target volume divided by the recycled volume per cube.
num_cubes_needed = volume_of_one_cube / total_recycled_volume_per_cube

# Print the final result, showing the numbers in the final equation
print("The calculation to find the number of cubes needed is:")
print("Total Volume of one Cube / Total Recycled Volume per Cube")
print("The final equation is:")
print(f"{int(volume_of_one_cube)} / {int(total_recycled_volume_per_cube)} = {int(num_cubes_needed)}")
<<<50>>>