import math

# Define the parameters based on the problem description
side_length = 10.0  # side length of the cube in mm
num_edges_chamfered = 4  # number of edges chamfered on one face
chamfer_hypotenuse = math.sqrt(2)  # width of the sloping edge in mm

# --- Step 1: Calculate the volume of one original cube ---
volume_of_one_cube = side_length ** 3

# --- Step 2: Calculate the volume of material removed from one edge ---
# The cross-section of the removed material is a right-angled isosceles triangle.
# Let 'a' be the length of the two equal sides.
# Using the Pythagorean theorem: a^2 + a^2 = chamfer_hypotenuse^2
# 2 * a^2 = (sqrt(2))^2 => 2 * a^2 = 2 => a^2 = 1 => a = 1 mm
a = math.sqrt(chamfer_hypotenuse**2 / 2.0)

# The volume of the removed triangular prism for one edge is (Area of triangle) * (length)
volume_removed_per_edge = (0.5 * a * a) * side_length

# --- Step 3: Calculate the total recycled volume per cube ---
total_volume_removed_per_cube = volume_removed_per_edge * num_edges_chamfered

# --- Step 4: Find the number of cubes needed ---
# N = Volume of one cube / Total volume removed per cube
# We use math.ceil because we need a whole number of cubes.
num_cubes_needed = math.ceil(volume_of_one_cube / total_volume_removed_per_cube)

# --- Step 5: Print the final equation and result ---
print("To find the number of chamfered cubes needed, we use the equation:")
print("Number of cubes = Volume of a new cube / Total recycled volume per cube")
print(f"Number of cubes = {volume_of_one_cube:.0f} / {total_volume_removed_per_cube:.0f}")
print(f"Number of cubes = {int(num_cubes_needed)}")

<<<50>>>