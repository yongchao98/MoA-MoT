import math
from fractions import Fraction

# Step 1: Define initial parameters of the cube and the chamfer process.
side_length = 10.0  # side length of the cube in mm
sloping_edge_width_squared = 2.0  # (sqrt(2))^2
num_chamfered_edges = 4

print("--- Calculation of Recycled Material and Number of Cubes ---")
print(f"\nA solid aluminum cube has a side length (s) of {side_length} mm.")

# Step 2: Calculate the volume of a single original cube.
volume_of_cube = side_length ** 3
print(f"The volume of one original cube is s^3 = {side_length}^3 = {volume_of_cube} mm^3.")

# Step 3: Determine the dimensions of the cut-off material.
# The chamfer cross-section is an isosceles right-angled triangle. The given sloping edge is the hypotenuse.
# Let the leg length be 'a'. Then a^2 + a^2 = (sloping_edge_width)^2.
# 2 * a^2 = sloping_edge_width_squared => a = sqrt(sloping_edge_width_squared / 2).
leg_a = math.sqrt(sloping_edge_width_squared / 2.0)
print(f"A chamfer is cut on {num_chamfered_edges} edges of one face.")
print(f"The leg 'a' of the triangular cut is {leg_a:.1f} mm deep into each face from the edge.")

# Step 4: Calculate the volume of the triangular prism removed for one edge.
triangle_area = 0.5 * leg_a ** 2
volume_per_edge = triangle_area * side_length
print(f"The volume of material removed for one edge (a triangular prism) is (1/2 * a^2) * s = (1/2 * {leg_a:.1f}^2) * {side_length} = {volume_per_edge:.1f} mm^3.")

# Step 5: Account for the overlaps at the corners.
# The volume of each corner overlap is a^3 / 3.
volume_per_overlap = (leg_a ** 3) / 3.0
total_volume_of_prisms = num_chamfered_edges * volume_per_edge
total_volume_of_overlaps = num_chamfered_edges * volume_per_overlap
net_recycled_volume = total_volume_of_prisms - total_volume_of_overlaps

# For exact fraction representation
net_recycled_volume_frac = Fraction(net_recycled_volume).limit_denominator()

print(f"\nTo find the total recycled volume, we must subtract the volume of the {num_chamfered_edges} corner overlaps.")
print(f"Volume of one corner overlap is a^3 / 3 = {leg_a:.1f}^3 / 3 = {volume_per_overlap:.4f} mm^3.")
print(f"The total recycled volume per cube is ({num_chamfered_edges} * {volume_per_edge:.1f}) - ({num_chamfered_edges} * {volume_per_overlap:.4f}) = {net_recycled_volume:.4f} mm^3 (exactly {net_recycled_volume_frac}).")

# Step 6: Calculate how many cubes are needed.
num_cubes_needed_float = volume_of_cube / net_recycled_volume
num_cubes_needed_integer = math.ceil(num_cubes_needed_float)

print(f"\nTo manufacture a new cube of {volume_of_cube} mm^3, we need to find N such that N * {net_recycled_volume_frac} >= {volume_of_cube}.")
print(f"The required number of cubes is {volume_of_cube} / {net_recycled_volume:.4f} = {num_cubes_needed_float:.4f}.")
print(f"Since the number of cubes must be a whole number, we round up.")
print(f"\n--- Final Answer ---")
print(f"The final number of chamfered cubes needed is: {num_cubes_needed_integer}")

print(f"\nThe final equation with all numbers is:")
print(f"N = ceil( (cube_side_length^3) / ({num_chamfered_edges} * (volume_per_edge) - {num_chamfered_edges} * (volume_per_overlap)) )")
print(f"N = ceil( ({side_length:.0f}^3) / ({num_chamfered_edges} * {(0.5 * leg_a**2 * side_length):.1f} - {num_chamfered_edges} * {(leg_a**3 / 3):.4f}) )")
print(f"N = ceil( {volume_of_cube:.0f} / ({total_volume_of_prisms:.1f} - {total_volume_of_overlaps:.4f}) )")
print(f"N = ceil( {volume_of_cube:.0f} / {net_recycled_volume:.4f} ) = {num_cubes_needed_integer}")
