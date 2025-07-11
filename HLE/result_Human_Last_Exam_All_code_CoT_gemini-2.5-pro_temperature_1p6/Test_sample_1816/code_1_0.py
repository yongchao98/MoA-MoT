import math

# Step 1: Define initial parameters
side_length = 10  # mm
num_chamfered_edges = 4
sloping_edge_width = math.sqrt(2)  # mm

# Step 2: Calculate the volume of a single cube
volume_of_cube = side_length ** 3

# Step 3: Calculate the volume of material removed from one chamfered edge
# The cross-section of the cut is an isosceles right-angled triangle.
# The sloping edge is the hypotenuse. Let 'a' be the length of the two equal sides.
# a^2 + a^2 = sloping_edge_width^2
# 2 * a^2 = (sqrt(2))^2 = 2
# a^2 = 1 => a = 1 mm
cut_depth = math.sqrt(sloping_edge_width**2 / 2)

# The shape of the removed material is a triangular prism.
# The area of its triangular base
area_of_triangular_cut = 0.5 * cut_depth * cut_depth

# The volume of the prism is its base area times its length (the cube's side length)
volume_per_chamfer = area_of_triangular_cut * side_length

# Step 4: Calculate the total recycled volume from one chamfered cube
total_recycled_volume_per_cube = volume_per_chamfer * num_chamfered_edges

# Step 5: Calculate the number of cubes needed to manufacture a new cube
# We need the accumulated recycled volume to be at least the volume of a new cube.
num_cubes_needed = math.ceil(volume_of_cube / total_recycled_volume_per_cube)

# --- Final Output ---
# Present the calculation step-by-step as an equation.
print("To find the number of cubes needed, we use the following equation:")
print("Number of Cubes = Volume of one Cube / Total Volume Recycled per Cube")
print("")
print("Volume of one Cube = (Side Length)^3")
print(f"                     = {side_length}^3 = {volume_of_cube} mm^3")
print("")
print("Total Volume Recycled per Cube = (Number of Edges) * (Volume per Chamfer)")
print(f"                               = {num_chamfered_edges} * {volume_per_chamfer} = {total_recycled_volume_per_cube} mm^3")
print("")
print("So, the final calculation is:")
print(f"Number of Cubes = {volume_of_cube} / {total_recycled_volume_per_cube}")
print(f"Number of Cubes = {int(num_cubes_needed)}")
<<<50>>>