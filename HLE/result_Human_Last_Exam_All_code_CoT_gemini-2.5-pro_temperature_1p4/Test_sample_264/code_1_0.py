import math

# Define the properties of the cubes
side_length = 1.0
shift_along_diagonal = 0.5

# Define the boundaries of the first cube (C1)
c1_min = 0.0
c1_max = c1_min + side_length

# The shift vector is applied to each coordinate
shift_vector_component = shift_along_diagonal

# Define the boundaries of the second, shifted cube (C2)
c2_min = c1_min + shift_vector_component
c2_max = c1_max + shift_vector_component

print("--- Step 1: Define the Cubes ---")
print(f"Cube 1 (C1) occupies the region where 0 <= x, y, z <= {c1_max:.1f}.")
print(f"Cube 2 (C2) is shifted by {shift_along_diagonal} along the main diagonal.")
print(f"C2 occupies the region where {c2_min:.1f} <= x, y, z <= {c2_max:.1f}.")

# Calculate the boundaries of the intersection
intersect_min = max(c1_min, c2_min)
intersect_max = min(c1_max, c2_max)

print("\n--- Step 2: Determine the Intersection ---")
print(f"The intersection is a smaller cube where the regions overlap.")
print(f"Its boundaries are from {intersect_min:.1f} to {intersect_max:.1f} on each axis.")

# Calculate the side length of the intersection cube
intersect_side_length = intersect_max - intersect_min
print(f"The side length of this intersection cube is {intersect_max:.1f} - {intersect_min:.1f} = {intersect_side_length:.1f}.")

# Calculate the volume of the intersection
volume = math.pow(intersect_side_length, 3)

print("\n--- Step 3: Calculate the Volume ---")
print("The volume of the intersection cube is its side length cubed.")
print(f"Volume = {intersect_side_length:.1f} * {intersect_side_length:.1f} * {intersect_side_length:.1f} = {volume}")