import math

# --- Parameters from the problem ---
side_length = 10  # in mm
# The width of the sloping edge (hypotenuse of the right-angled triangle cut)
chamfer_hypotenuse = math.sqrt(2)
num_chamfered_edges = 4

# --- Step 1: Calculate the volume of a standard cube ---
volume_cube = side_length ** 3

# --- Step 2: Calculate the volume of one chamfered-off piece ---
# The cross-section of the chamfered piece is a right-angled isosceles triangle.
# Let 'x' be the length of the two equal sides.
# By Pythagorean theorem: x^2 + x^2 = chamfer_hypotenuse^2
# 2 * x^2 = 2, so x^2 = 1, and x = 1.
# This 'x' is the cut depth along each face.
cut_depth = 1 # in mm

# The volume of the cut-off piece (a triangular prism) is the area of the triangle * length.
# Area of the right-angled triangle cross-section
triangle_area = 0.5 * cut_depth * cut_depth
# Volume of the material removed from one edge
volume_per_chamfer = triangle_area * side_length

# --- Step 3: Calculate the total recycled material from one cube ---
total_recycled_volume_per_cube = volume_per_chamfer * num_chamfered_edges

# --- Step 4: Calculate the number of cubes needed ---
# This is the total volume needed divided by the volume recycled from each cube.
num_cubes_needed = volume_cube / total_recycled_volume_per_cube

# --- Final Output ---
print("To find the number of chamfered cubes needed, we use the following equation:")
print("Number of cubes = (Total Volume of one Cube) / (Recycled Volume per Cube)")
print("The final calculation is:")
print(f"{int(volume_cube)} / {int(total_recycled_volume_per_cube)} = {int(num_cubes_needed)}")