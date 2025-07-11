# Define the lengths of the initial structure (Iteration 0)
l0_trunk = 40
l0_branch = 20

# The scaling factor is the ratio of the initial branch length to the initial trunk length.
scaling_factor = l0_branch / l0_trunk # 20 / 40 = 0.5

# Calculate the lengths of the segments for the structure added in Iteration 1
# This structure replaces a branch of length 20.
l1_trunk = scaling_factor * l0_trunk
l1_branch = scaling_factor * l0_branch

# Calculate the lengths of the segments for the structure added in Iteration 2
# This structure replaces a branch of length l1_branch (10).
# Its trunk is scaled from the trunk of the previous structure (l1_trunk).
l2_trunk = scaling_factor * l1_trunk
l2_branch = scaling_factor * l1_branch

# The white path consists of the main trunk, the trunk from the 1st iteration,
# the trunk from the 2nd iteration, and the two final branches from the 2nd iteration.
segment1 = l0_trunk
segment2 = l1_trunk
segment3 = l2_trunk
segment4 = l2_branch
segment5 = l2_branch

# Calculate the total length
total_length = segment1 + segment2 + segment3 + segment4 + segment5

# Print the final equation with each number
print(f"The total length is the sum of the following segments:")
print(f"{int(segment1)} (main trunk) + {int(segment2)} (iter 1 trunk) + {int(segment3)} (iter 2 trunk) + {int(segment4)} (iter 2 branch) + {int(segment5)} (iter 2 branch)")

# Print the final answer
print(f"Total length = {int(total_length)}")
