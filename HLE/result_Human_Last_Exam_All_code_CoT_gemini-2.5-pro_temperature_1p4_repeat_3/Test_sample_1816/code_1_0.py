import math

# --- Step 1: Define initial parameters ---
cube_side_length = 10  # in mm
chamfer_slope_width = math.sqrt(2)  # in mm
num_chamfered_edges = 4

# --- Step 2: Calculate the volume of a new cube ---
volume_of_new_cube = cube_side_length ** 3

# --- Step 3: Calculate the volume of material removed from one cube ---
# The cross-section of the chamfer is a right-angled isosceles triangle.
# The chamfer_slope_width is the hypotenuse. Let the other two sides be 'a'.
# a^2 + a^2 = (sqrt(2))^2  => 2a^2 = 2 => a = 1 mm.
cut_depth = 1 # mm
# The area of the triangular cross-section
chamfer_cross_section_area = 0.5 * cut_depth * cut_depth
# The volume of one chamfer (a triangular prism)
volume_per_chamfer = chamfer_cross_section_area * cube_side_length
# Total recycled volume from one cube
recycled_volume_per_cube = volume_per_chamfer * num_chamfered_edges

# --- Step 4: Calculate the number of cubes needed ---
# We need enough recycled material to equal the volume of one new cube.
num_cubes_needed = volume_of_new_cube / recycled_volume_per_cube

# --- Step 5: Print the results and the final equation ---
print(f"The volume of one new cube is: {volume_of_new_cube} mm^3.")
print(f"The total volume of recycled material from one chamfered cube is: {recycled_volume_per_cube} mm^3.")
print("\nTo find the number of cubes needed, we use the equation:")
print("Number of Cubes = (Volume of a New Cube) / (Recycled Volume per Cube)")

# Print the final equation with the calculated numbers
print(f"\n{int(volume_of_new_cube)} / {int(recycled_volume_per_cube)} = {int(num_cubes_needed)}")

print(f"\nTherefore, {int(num_cubes_needed)} chamfered cubes are needed to manufacture one new cube.")

<<<50>>>