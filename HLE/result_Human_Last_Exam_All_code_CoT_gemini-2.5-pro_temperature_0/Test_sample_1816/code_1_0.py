import math

# Step 1: Define the properties of the cube and the chamfer process.
side_length = 10  # in mm
num_chamfered_edges = 4
chamfer_hypotenuse = math.sqrt(2) # in mm

# Step 2: Calculate the volume of a single, original cube.
volume_cube = side_length ** 3

# Step 3: Calculate the volume of material removed by one chamfer.
# The cross-section of the chamfer cut is a right-angled isosceles triangle.
# Let 'a' be the length of the two equal sides.
# According to the Pythagorean theorem: a^2 + a^2 = chamfer_hypotenuse^2
# 2 * a^2 = (sqrt(2))^2 = 2  => a^2 = 1 => a = 1 mm.
cut_side_length = 1.0
# The area of this triangular cross-section is (1/2) * base * height.
cut_cross_section_area = 0.5 * cut_side_length * cut_side_length
# The volume of the removed prism is the cross-section area times the length of the edge.
volume_per_chamfer = cut_cross_section_area * side_length

# Step 4: Calculate the total volume of recycled material from one cube.
total_recycled_volume_per_cube = num_chamfered_edges * volume_per_chamfer

# Step 5: Calculate how many cubes are needed to make one new cube.
# This is the volume of a full cube divided by the recycled volume per cube.
num_cubes_needed = volume_cube / total_recycled_volume_per_cube

# Step 6: Print the final equation and the result.
print("To find the number of cubes needed, we set up the following equation:")
print("Number of Cubes * (Recycled Volume per Cube) = Volume of one New Cube")
print(f"Number of Cubes * {total_recycled_volume_per_cube:.0f} = {volume_cube:.0f}")
print("Solving for the number of cubes:")
print(f"Number of Cubes = {volume_cube:.0f} / {total_recycled_volume_per_cube:.0f} = {num_cubes_needed:.0f}")

print("\nTherefore, 50 chamfered cubes are needed to recycle enough material for one new cube.")
<<<50>>>