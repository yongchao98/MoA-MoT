import math

# --- 1. Define initial parameters ---
side_length = 10  # in mm
chamfer_slope_width = math.sqrt(2) # in mm
num_chamfered_edges = 4

print(f"Cube side length: {side_length} mm")
print(f"Number of chamfered edges: {num_chamfered_edges}")
print("-" * 30)

# --- 2. Calculate the volume of a single cube ---
cube_volume = side_length ** 3
print(f"Volume of one complete cube: {side_length} * {side_length} * {side_length} = {cube_volume} mm^3")
print("-" * 30)

# --- 3. Calculate the volume of removed material per cube ---
# For a 45-degree right triangle, a^2 + a^2 = c^2, so 2a^2 = c^2
# a = sqrt(c^2 / 2)
cut_depth_a = math.sqrt(chamfer_slope_width**2 / 2)

# Volume of one removed prism = Area of triangular base * length
triangle_area = 0.5 * cut_depth_a * cut_depth_a
removed_volume_per_edge = triangle_area * side_length
print(f"Cut depth from each face of the edge: {cut_depth_a:.1f} mm")
print(f"Volume removed from a single edge: {triangle_area:.1f} * {side_length} = {removed_volume_per_edge:.1f} mm^3")

# Total volume removed from one cube
total_removed_volume = removed_volume_per_edge * num_chamfered_edges
print(f"Total recycled volume from one cube: {num_chamfered_edges} * {removed_volume_per_edge:.1f} = {total_removed_volume:.1f} mm^3")
print("-" * 30)

# --- 4. Calculate how many cubes are needed ---
# This must be a whole number, so we use ceiling division or simple division if it's exact.
required_cubes = cube_volume / total_removed_volume

print("To find the number of chamfered cubes needed, we use the following equation:")
print("Required Cubes = (Volume of one Cube) / (Total Recycled Volume from one Cube)")
print("\nFinal Calculation:")
# Final print with each number from the equation
print(f"{cube_volume:.0f} / {total_removed_volume:.0f} = {required_cubes:.0f}")

print(f"\nTherefore, {int(required_cubes)} chamfered cubes are needed to manufacture another cube from the recycled material.")
