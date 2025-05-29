def transform_grid(input_grid):
    # Find the first non-zero element
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    if first_non_zero_index is None:
        return input_grid  # No non-zero elements, return as is

    # Find the last non-zero element
    last_non_zero_index = len(input_grid) - next((i for i, x in enumerate(reversed(input_grid)) if x != 0), None) - 1

    # Extract the non-zero sequence
    non_zero_sequence = input_grid[first_non_zero_index:last_non_zero_index + 1]

    # Calculate the number of leading zeros
    leading_zeros = first_non_zero_index

    # Create the new output grid
    output_grid = [0] * leading_zeros + [0] * leading_zeros + non_zero_sequence + [0] * (len(input_grid) - len(non_zero_sequence) - 2 * leading_zeros)

    return output_grid

# Test input
test_input = [0, 0, 8, 8, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
output = transform_grid(test_input)
print(output)