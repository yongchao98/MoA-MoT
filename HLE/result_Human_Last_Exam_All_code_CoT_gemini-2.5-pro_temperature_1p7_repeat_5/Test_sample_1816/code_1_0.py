import math

# Step 1: Define initial parameters
side_length = 10  # in mm
chamfer_width = math.sqrt(2)  # in mm
num_chamfered_edges = 4

print(f"A solid cube has a side length of {side_length} mm.")

# Step 2: Calculate the volume of the original cube
volume_cube = side_length ** 3
print(f"The volume of one cube is {side_length}^3 = {volume_cube} mm^3.")
print("-" * 20)

# Step 3: Calculate the volume of material removed by one chamfer
# The cross-section is a right-angled isosceles triangle with hypotenuse = chamfer_width.
# Let the other two equal sides be 'x'. Then x^2 + x^2 = chamfer_width^2
# 2*x^2 = (sqrt(2))^2 => 2*x^2 = 2 => x = 1 mm
cut_depth = 1
# Area of the triangular cross-section
area_cut_cross_section = 0.5 * cut_depth * cut_depth
# Volume of one chamfer cut (a triangular prism)
volume_one_chamfer = area_cut_cross_section * side_length
print(f"The volume of material removed from one edge (one chamfer) is {volume_one_chamfer} mm^3.")

# Step 4: Calculate the total recycled volume per cube
total_recycled_volume_per_cube = volume_one_chamfer * num_chamfered_edges
print(f"Since {num_chamfered_edges} edges are chamfered, the total recycled material from one cube is {num_chamfered_edges} * {volume_one_chamfer} = {total_recycled_volume_per_cube} mm^3.")
print("-" * 20)

# Step 5: Calculate how many chamfered cubes are needed to make a new one
# Use math.ceil to ensure we have enough material (though in this case it's a perfect integer)
num_cubes_needed = math.ceil(volume_cube / total_recycled_volume_per_cube)

# Final Output
print("To find the number of chamfered cubes needed to make a new cube, we divide the volume of a full cube by the recycled volume per cube.")
print(f"Final Calculation: {int(volume_cube)} / {int(total_recycled_volume_per_cube)} = {int(num_cubes_needed)}")
print("\nTherefore, 50 chamfered cubes are needed to make the accumulated recycled material enough for manufacturing another cube.")

print(f"\n<<<{int(num_cubes_needed)}>>>")