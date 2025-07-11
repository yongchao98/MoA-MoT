import math

# Define the initial parameters of the cube and the chamfer process.
side_length = 10  # in mm
chamfer_hypotenuse = math.sqrt(2)  # width of the sloping edge in mm
num_chamfered_edges = 4

# --- Calculation Steps ---

# Step 1: Calculate the volume of a full, original cube.
volume_full_cube = side_length**3

# Step 2: Determine the dimensions of the material removed by one chamfer.
# The chamfer cut removes a triangular prism. The base of this prism is a
# right-angled isosceles triangle. The given chamfer width is the hypotenuse.
# For a right-angled isosceles triangle with equal sides 'a': a^2 + a^2 = hypotenuse^2
# 2 * a^2 = (sqrt(2))^2 => 2 * a^2 = 2 => a^2 = 1 => a = 1 mm.
cut_depth = math.sqrt(chamfer_hypotenuse**2 / 2)

# Step 3: Calculate the volume of one removed triangular prism (one chamfer).
# Area of the triangular base = 0.5 * base * height = 0.5 * cut_depth * cut_depth
triangle_base_area = 0.5 * cut_depth * cut_depth
# Volume of one chamfer = Area of base * length of the prism (which is the side_length)
volume_one_chamfer = triangle_base_area * side_length

# Step 4: Calculate the total volume of recycled material from one cube.
total_recycled_volume_per_cube = volume_one_chamfer * num_chamfered_edges

# Step 5: Calculate how many chamfered cubes are needed.
# This is the volume of a full cube divided by the recycled volume per cube.
num_cubes_needed = volume_full_cube / total_recycled_volume_per_cube

# --- Output the result ---

print(f"The volume of one new cube is {volume_full_cube} mm^3.")
print(f"The volume of recycled material from one chamfered cube is {total_recycled_volume_per_cube} mm^3.")
print("\nTo find the number of cubes needed, the following equation is solved:")
print(f"Number of Cubes = (Volume of a New Cube) / (Recycled Volume per Cube)")
print(f"Number of Cubes = {volume_full_cube:.0f} / {total_recycled_volume_per_cube:.0f} = {num_cubes_needed:.0f}")
<<<50>>>