import math

# --- Define initial parameters ---
side_length = 10  # in mm
chamfer_width = math.sqrt(2) # in mm, this is the hypotenuse
num_chamfered_edges = 4

# --- Step 1: Calculate the volume of the original cube ---
volume_cube = side_length**3

# --- Step 2: Calculate the volume of a single chamfer ---
# The cross-section is a right-angled isosceles triangle.
# Let 'a' be the side length of the triangle legs on the cube's faces.
# a^2 + a^2 = chamfer_width^2  => 2*a^2 = 2 => a = 1 mm.
chamfer_leg = math.sqrt(chamfer_width**2 / 2)

# The volume of one removed chamfer (a triangular prism) is
# (cross-sectional area) * (length of the edge).
volume_one_chamfer = (0.5 * chamfer_leg * chamfer_leg) * side_length

# --- Step 3: Calculate the total recycled volume per cube ---
total_volume_removed = volume_one_chamfer * num_chamfered_edges

# --- Step 4: Calculate how many cubes are needed ---
# This is the volume of a full cube divided by the recycled volume per cube.
num_cubes_needed = volume_cube / total_volume_removed

# --- Print the final equation with all the numbers ---
print("The final calculation is based on the equation:")
print("Number of Cubes = (Volume of one Cube) / (Volume of one Chamfer * Number of Chamfers)")
print("\nPlugging in the numbers:")
# The final line prints the equation with the calculated values
print(f"{int(num_cubes_needed)} = {int(volume_cube)} / ({int(volume_one_chamfer)} * {int(num_chamfered_edges)})")