import math

# Define the initial parameters of the problem
side_length = 10.0  # mm
num_chamfered_edges = 4
chamfer_hypotenuse = math.sqrt(2.0)

# Step 1: Calculate the volume of a full cube
volume_of_full_cube = side_length ** 3

# Step 2: Calculate the volume of material from one chamfer cut
# The cross-section is a right-angled isosceles triangle.
# Let 'a' be the length of the two equal sides.
# a^2 + a^2 = chamfer_hypotenuse^2  => 2*a^2 = 2 => a^2 = 1 => a = 1
cut_side_length = 1.0

# Area of the triangular cross-section = 0.5 * base * height
triangle_area = 0.5 * cut_side_length * cut_side_length

# Volume of one chamfer cut (triangular prism) = area * length
volume_of_one_chamfer = triangle_area * side_length

# Step 3: Calculate the total recycled volume from one chamfered cube
total_recycled_volume = num_chamfered_edges * volume_of_one_chamfer

# Step 4: Calculate how many cubes are needed
num_cubes_needed = volume_of_full_cube / total_recycled_volume

print("This script calculates the number of chamfered cubes required to recycle enough material for one new cube.")
print("Here is the breakdown of the calculation:")

# Present the final equation with all the numbers
print(f"\nNumber of Cubes = (Volume of Full Cube) / (Total Recycled Volume from one Cube)")
print(f"= ({side_length:.0f} * {side_length:.0f} * {side_length:.0f}) / ({num_chamfered_edges:.0f} * (0.5 * {cut_side_length:.0f} * {cut_side_length:.0f} * {side_length:.0f}))")
print(f"= {volume_of_full_cube:.0f} / ({num_chamfered_edges:.0f} * {volume_of_one_chamfer:.0f})")
print(f"= {volume_of_full_cube:.0f} / {total_recycled_volume:.0f}")
print(f"= {num_cubes_needed:.0f}")
print(f"\nTherefore, {int(num_cubes_needed)} chamfered cubes are needed to manufacture another cube.")

print("\n<<<50>>>")