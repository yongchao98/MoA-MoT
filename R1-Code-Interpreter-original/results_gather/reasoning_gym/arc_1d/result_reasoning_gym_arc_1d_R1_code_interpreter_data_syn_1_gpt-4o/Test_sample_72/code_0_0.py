def transform_grid(input_grid):
    # Find the first non-zero element
    start_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    if start_index is None:
        return input_grid  # If all zeros, return as is

    # Find the end of the non-zero sequence
    end_index = next((i for i, x in enumerate(input_grid[start_index:], start=start_index) if x == 0), len(input_grid))

    # Extract the non-zero sequence
    non_zero_sequence = input_grid[start_index:end_index]

    # Create the output grid with zeros
    output_grid = [0] * len(input_grid)

    # Place the non-zero sequence shifted by 3 positions
    output_start_index = start_index + 3
    output_grid[output_start_index:output_start_index + len(non_zero_sequence)] = non_zero_sequence

    return output_grid

# Test input
test_input = [0, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
output = transform_grid(test_input)
print(output)