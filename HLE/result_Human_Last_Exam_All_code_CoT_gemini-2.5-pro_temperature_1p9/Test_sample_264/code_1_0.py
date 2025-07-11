# Define the properties of the cubes and the shift
cube_side_length = 1.0
shift = 0.5

# Let's define the first cube (C1) to be in the region [0, 1] on each axis.
c1_min = 0.0
c1_max = c1_min + cube_side_length

# The second cube (C2) is shifted by 'shift' along each axis.
c2_min = c1_min + shift
c2_max = c1_max + shift

# The intersection is a smaller cube. Its boundaries are the maximum of the minimums
# and the minimum of the maximums.
intersection_min = max(c1_min, c2_min)
intersection_max = min(c1_max, c2_max)

# Calculate the side length of the intersection cube.
intersection_side_length = intersection_max - intersection_min

# Calculate the volume of the intersection.
volume = intersection_side_length ** 3

# Print the final equation as requested.
print(f"The volume of the intersection is calculated as the cube of its side length.")
print(f"Side length of intersection = {intersection_side_length}")
print(f"Volume = {intersection_side_length} * {intersection_side_length} * {intersection_side_length} = {volume}")