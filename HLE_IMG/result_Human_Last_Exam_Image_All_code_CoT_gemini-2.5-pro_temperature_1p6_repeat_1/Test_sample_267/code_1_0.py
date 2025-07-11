# Initial lengths of the base structure (iteration 0)
l0_trunk = 40
l0_branch = 20

# The white path traces the main trunk and the 'newest' path of growth.

# The first segment of the white path is the main trunk.
segment1_length = l0_trunk

# In the 1st iteration, a branch of length 20 is replaced by a new 'Y' structure.
# The trunk of this new structure has length 20. This is the second segment of the white path.
segment2_length = l0_branch
# The branches of this structure would be 20 * (20/40) = 10.

# In the 2nd iteration, a branch of length 10 is replaced by another 'Y' structure.
# The trunk of this newest structure has length 10. This is the third segment.
segment3_length = segment2_length * (l0_branch / l0_trunk)

# The final two parts of the white path are the two branches from the 2nd iteration.
# Their length is their trunk's length (10) times the scaling ratio.
segment4_and_5_length = segment3_length * (l0_branch / l0_trunk)

# Calculate the total length
total_length = segment1_length + segment2_length + segment3_length + segment4_and_5_length + segment4_and_5_length

# Print the equation with the calculated values
print(f"The total length of the white path is the sum of its segments:")
print(f"{int(segment1_length)} (main trunk) + {int(segment2_length)} (iter 1 trunk) + {int(segment3_length)} (iter 2 trunk) + {int(segment4_and_5_length)} (iter 2 branch) + {int(segment4_and_5_length)} (iter 2 branch) = {int(total_length)}")