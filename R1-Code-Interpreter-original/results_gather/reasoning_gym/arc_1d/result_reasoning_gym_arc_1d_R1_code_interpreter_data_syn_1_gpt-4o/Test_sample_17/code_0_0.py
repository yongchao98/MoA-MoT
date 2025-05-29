def transform_grid(input_grid):
    # Find the first non-zero element
    start_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    if start_index is None:
        return input_grid  # No non-zero element found, return as is

    # Find the last non-zero element
    end_index = next((i for i, x in enumerate(reversed(input_grid)) if x != 0), None)
    if end_index is None:
        return input_grid  # No non-zero element found, return as is

    end_index = len(input_grid) - end_index - 1

    # Extract the sequence of non-zero elements
    non_zero_sequence = input_grid[start_index:end_index + 1]

    # Shift the sequence one position to the left
    shifted_sequence = non_zero_sequence[1:] + [0]

    # Construct the new grid
    new_grid = input_grid[:start_index] + shifted_sequence + input_grid[end_index + 1:]

    return new_grid

# Test input
test_input = [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0]
output_grid = transform_grid(test_input)
print(output_grid)