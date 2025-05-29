def transform_grid(input_grid):
    # Find the first non-zero index
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), len(input_grid))
    # Extract the non-zero sequence
    non_zero_sequence = input_grid[first_non_zero_index:]
    # Calculate the number of zeros to append
    num_zeros = len(input_grid) - len(non_zero_sequence)
    # Create the output grid
    output_grid = non_zero_sequence + [0] * num_zeros
    return output_grid

# Test input
test_input = [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]
# Transform the test input
output = transform_grid(test_input)
print(output)