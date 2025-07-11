import math

# Step 1: Define initial parameters and calculate the volume of the original cube.
side_length = 10  # in mm
volume_of_cube = side_length ** 3

# Step 2: Calculate the volume of material removed by a single chamfer.
# The chamfer is a 45-degree cut, creating a right-angled isosceles triangle cross-section.
# The sloping edge width is the hypotenuse of this triangle.
hypotenuse = math.sqrt(2)  # in mm

# Using Pythagorean theorem (a^2 + a^2 = c^2) to find the leg 'a' of the triangle.
# 2 * a^2 = hypotenuse^2
# a = sqrt(hypotenuse^2 / 2)
leg_a_squared = (hypotenuse**2) / 2
leg_a = math.sqrt(leg_a_squared)

# The area of the triangular cross-section is (1/2) * base * height.
area_of_triangle = 0.5 * leg_a * leg_a

# The volume of the removed prism (one chamfer) is the area * length.
volume_per_chamfer = area_of_triangle * side_length

# Step 3: Calculate the total recycled volume from one cube.
# Four edges on one face are chamfered.
num_chamfers = 4
total_volume_removed_per_cube = num_chamfers * volume_per_chamfer

# Step 4: Calculate the number of cubes needed to make one new cube.
# This is the volume of a full cube divided by the recycled volume per cube.
# We use math.ceil to ensure we have enough material.
num_cubes_needed = math.ceil(volume_of_cube / total_volume_removed_per_cube)

# Output the final equation with all the calculated numbers.
print("To find the number of chamfered cubes needed, we use the following equation:")
print("Number of Cubes = Total Volume of One Cube / (Number of Chamfers * Volume per Chamfer)")
print(f"Number of Cubes = {volume_of_cube} mm^3 / ({num_chamfers} * {volume_per_chamfer} mm^3)")
print(f"Number of Cubes = {volume_of_cube} / {total_volume_removed_per_cube}")
print(f"Result: {int(num_cubes_needed)}")

<<<50>>>