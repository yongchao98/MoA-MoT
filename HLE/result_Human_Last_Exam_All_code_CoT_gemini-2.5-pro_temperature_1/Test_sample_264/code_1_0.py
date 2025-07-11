# Define the properties of the cubes based on the problem description.
cube_side_length = 1.0
shift_distance = 0.5

# The first cube occupies the interval [0, 1] on each axis.
# The second cube is shifted, occupying the interval [0.5, 1.5] on each axis.

# To find the intersection, we determine the overlapping interval for any given axis.
# The start of the overlap is the maximum of the two starting points.
overlap_start = max(0.0, shift_distance)

# The end of the overlap is the minimum of the two ending points.
overlap_end = min(cube_side_length, cube_side_length + shift_distance)

# The side length of the resulting intersection cube is the length of this overlap.
intersection_side_length = overlap_end - overlap_start

# The volume of the intersection is its side length cubed.
volume = intersection_side_length ** 3

print("Step 1: The two cubes are defined as:")
print(f"Cube 1: 0 <= x,y,z <= {cube_side_length}")
print(f"Cube 2: {shift_distance} <= x,y,z <= {cube_side_length + shift_distance}")
print("-" * 20)
print("Step 2: The intersection is a smaller cube with side length:")
print(f"Side Length = min({cube_side_length}, {cube_side_length + shift_distance}) - max(0.0, {shift_distance}) = {intersection_side_length}")
print("-" * 20)
print("Step 3: The volume is the side length cubed. The final equation is:")
# The final equation with each number printed, as requested.
print(f"{intersection_side_length} * {intersection_side_length} * {intersection_side_length} = {volume}")