import math

# --- Problem Parameters ---
side_length = 10  # in mm
chamfer_hypotenuse = math.sqrt(2) # in mm
num_chamfered_edges = 4

# Step 1: Calculate the volume of the original cube
volume_cube = side_length**3
print(f"First, we calculate the volume of a single aluminum cube.")
print(f"Volume of one cube = {side_length} * {side_length} * {side_length} = {int(volume_cube)} mm^3\n")

# Step 2: Calculate the volume of material removed from one chamfered edge.
# The cross-section of the removed material is an isosceles right-angled triangle.
# The hypotenuse is sqrt(2) mm. Let the other two sides be 'a'.
# a^2 + a^2 = (sqrt(2))^2  =>  2a^2 = 2  =>  a = 1 mm.
chamfer_leg_side = 1 # in mm

# The area of this triangular cross-section is (1/2) * base * height.
area_triangle = 0.5 * chamfer_leg_side * chamfer_leg_side

# The volume of the removed prism is the cross-section area times its length (the cube's side length).
volume_one_chamfer = area_triangle * side_length

print(f"Next, we calculate the volume of the material cut off from one edge.")
print(f"The cross-section is a right-angled triangle with two equal sides of {chamfer_leg_side} mm.")
print(f"Volume of one chamfer cut = (1/2 * {chamfer_leg_side} * {chamfer_leg_side}) * {side_length} = {int(volume_one_chamfer)} mm^3\n")

# Step 3: Calculate the total recycled volume from one cube.
total_recycled_volume = volume_one_chamfer * num_chamfered_edges
print(f"One cube has {num_chamfered_edges} edges chamfered.")
print(f"Total recycled volume per cube = {num_chamfered_edges} * {int(volume_one_chamfer)} = {int(total_recycled_volume)} mm^3\n")

# Step 4: Calculate how many cubes are needed.
# This is the volume of a new cube divided by the recycled volume per cube.
# We use math.ceil to ensure we have enough material, in case the division results in a fraction.
num_cubes_needed = math.ceil(volume_cube / total_recycled_volume)

print(f"Finally, we find the number of chamfered cubes needed to make one new cube.")
print(f"Number of cubes needed = (Volume of one cube) / (Recycled volume per cube)")
print(f"Equation: {int(volume_cube)} / {int(total_recycled_volume)} = {int(num_cubes_needed)}")
print(f"\nTherefore, {int(num_cubes_needed)} chamfered cubes are needed to recycle enough material to make one new cube.")

<<<50>>>