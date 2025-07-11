# Step 1 & 2: Define initial parameters and scaling rule.
# The initial structure has a main trunk and two branches.
initial_trunk_length = 40
initial_branch_length = 20

# The scaling rule is that a branch of length L is replaced by a new "Y" structure.
# The trunk of the new "Y" has length L.
# The branches of the new "Y" have length L * (initial_branch_length / initial_trunk_length), which is L * 0.5.

# Step 3 & 4: Trace the iterations and identify the white branches.
# The white path in the final image consists of several segments.

# Segment 1: The main trunk from the original structure.
white_segment_1 = initial_trunk_length

# Segment 2: The trunk of the "Y" created in the 1st iteration.
# This replaced a branch of length 20. So, its trunk's length is 20.
white_segment_2 = initial_branch_length

# Segment 3: The trunk of the "Y" created in the 2nd iteration.
# This replaced a branch from the 1st iteration's "Y". Those branches had length 20 * 0.5 = 10.
# So, the trunk of this newest "Y" has length 10.
white_segment_3_trunk = initial_branch_length * 0.5

# Segments 4 & 5: The two branches of the "Y" created in the 2nd iteration.
# Their length is half of their trunk's length (10 * 0.5 = 5).
white_segment_3_branch = white_segment_3_trunk * 0.5

# Step 5: Calculate the total length.
total_length = white_segment_1 + white_segment_2 + white_segment_3_trunk + white_segment_3_branch + white_segment_3_branch

# Print the final calculation, showing each component.
print(f"The total length of the white branches is the sum of its parts:")
print(f"Main trunk: {int(white_segment_1)}")
print(f"Trunk from 1st iteration: {int(white_segment_2)}")
print(f"Trunk from 2nd iteration: {int(white_segment_3_trunk)}")
print(f"Branches from 2nd iteration: {int(white_segment_3_branch)} and {int(white_segment_3_branch)}")
print(f"Final Equation: {int(white_segment_1)} + {int(white_segment_2)} + {int(white_segment_3_trunk)} + {int(white_segment_3_branch)} + {int(white_segment_3_branch)} = {int(total_length)}")
