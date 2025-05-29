# Define the input grid
input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 8, 9, 2, 8, 8, 6, 3, 9, 8, 2, 6, 4, 7, 7]

# Find the first non-zero index
first_non_zero_index = next(i for i, x in enumerate(input_grid) if x != 0)

# Extract the non-zero sequence
non_zero_sequence = input_grid[first_non_zero_index:]

# Shift the sequence to the left and move the first element to the end
shifted_sequence = non_zero_sequence[1:] + [non_zero_sequence[0]]

# Construct the output grid
output_grid = shifted_sequence + [0] * (len(input_grid) - len(shifted_sequence))

# Print the output grid
print(output_grid)