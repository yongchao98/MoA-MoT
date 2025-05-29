def transform_grid(input_grid):
    # Find the first non-zero sequence
    first_non_zero = next((i for i, x in enumerate(input_grid) if x != 0), None)
    if first_non_zero is None:
        return input_grid  # No non-zero element found, return as is

    # Find the end of the first non-zero sequence
    end_non_zero = next((i for i, x in enumerate(input_grid[first_non_zero:], start=first_non_zero) if x == 0), len(input_grid))

    # Extend the non-zero sequence by two
    non_zero_sequence = input_grid[first_non_zero:end_non_zero]
    extended_sequence = non_zero_sequence + non_zero_sequence[:2]

    # Shift the zero sequence by two positions
    zero_sequence = input_grid[end_non_zero:]
    shifted_zero_sequence = [0, 0] + zero_sequence[:-2]

    # Combine the sequences
    output_grid = extended_sequence + shifted_zero_sequence
    return output_grid

# Test input
input_grid = [4, 4, 4, 4, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4]
output_grid = transform_grid(input_grid)
print(output_grid)