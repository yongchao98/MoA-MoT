def transform_grid(input_grid):
    # Find the first non-zero element
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    if first_non_zero_index is None:
        return input_grid  # No non-zero elements, return as is

    # Find the last non-zero element
    last_non_zero_index = len(input_grid) - 1 - next((i for i, x in enumerate(reversed(input_grid)) if x != 0), None)

    # Determine the shift amount
    shift_amount = 4

    # Create a new grid with zeros
    output_grid = [0] * len(input_grid)

    # Shift the sequence to the right
    for i in range(first_non_zero_index, last_non_zero_index + 1):
        if input_grid[i] != 0:
            new_position = i + shift_amount
            if new_position < len(input_grid):
                output_grid[new_position] = input_grid[i]

    return output_grid

# Test input
test_input = [9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9]
output = transform_grid(test_input)
print(output)