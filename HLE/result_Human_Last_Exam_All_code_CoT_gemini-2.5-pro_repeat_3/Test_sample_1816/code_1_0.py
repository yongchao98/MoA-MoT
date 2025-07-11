import math

# Step 1: Define parameters and calculate the volume of the original cube.
side_length = 10.0
volume_cube = side_length ** 3

print(f"The volume of one original cube is {side_length} * {side_length} * {side_length} = {volume_cube} mm^3.")
print("-" * 20)

# Step 2: Calculate the volume of a single chamfer cut.
# The chamfer has a sloping edge (hypotenuse) of sqrt(2) and is at 45 degrees.
# The cross-section is an isosceles right-angled triangle.
# Let 'a' be the two equal sides of the triangle. a^2 + a^2 = (sqrt(2))^2 => 2a^2 = 2 => a = 1.
cut_depth_a = 1.0
cross_sectional_area = 0.5 * cut_depth_a * cut_depth_a
volume_one_chamfer = cross_sectional_area * side_length

print(f"The volume of material removed by one chamfer is {cross_sectional_area} * {side_length} = {volume_one_chamfer} mm^3.")

# Step 3: Calculate the total recycled volume from one cube.
num_chamfered_edges = 4
total_recycled_volume = volume_one_chamfer * num_chamfered_edges

print(f"The total recycled volume from one cube (with {num_chamfered_edges} chamfers) is {volume_one_chamfer} * {num_chamfered_edges} = {total_recycled_volume} mm^3.")
print("-" * 20)

# Step 4: Determine the number of cubes needed.
num_cubes_needed = volume_cube / total_recycled_volume

print("To find how many chamfered cubes are needed to make a new cube, we form the equation:")
print(f"Number of Cubes = (Volume of a Cube) / (Recycled Volume per Cube)")
print(f"Number of Cubes = {volume_cube} / {total_recycled_volume} = {int(num_cubes_needed)}")
<<<50>>>