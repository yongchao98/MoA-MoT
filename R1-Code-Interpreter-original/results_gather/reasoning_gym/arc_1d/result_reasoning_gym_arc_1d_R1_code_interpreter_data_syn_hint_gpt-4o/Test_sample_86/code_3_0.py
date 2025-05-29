# Test input grid
input_grid = [8, 1, 2, 0, 0, 0, 0, 0, 0, 0, 8, 7, 2, 8, 4, 5, 1, 8, 2]

# Identify the non-zero sequence
non_zero_sequence = [x for x in input_grid if x != 0]

# Determine the position to shift the sequence
# Since the sequence is at the start, move it to the end
output_grid = [0] * (len(input_grid) - len(non_zero_sequence)) + non_zero_sequence

# Print the output grid
print(output_grid)