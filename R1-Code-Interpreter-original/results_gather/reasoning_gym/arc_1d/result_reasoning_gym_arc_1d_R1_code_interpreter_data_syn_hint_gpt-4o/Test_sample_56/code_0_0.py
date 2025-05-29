def transform_grid(input_grid):
    # Find the non-zero sequence
    non_zero_sequence = [x for x in input_grid if x != 0]
    # Create an output grid of zeros
    output_grid = [0] * len(input_grid)
    # Determine the starting position for the shifted sequence
    start_position = len(input_grid) - len(non_zero_sequence) - 7
    # Place the non-zero sequence in the output grid
    for i, value in enumerate(non_zero_sequence):
        output_grid[start_position + i] = value
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0]
# Transform the grid
output_grid = transform_grid(input_grid)
print(output_grid)