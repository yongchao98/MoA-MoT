def transform_grid(input_grid):
    output_grid = [row[:] for row in input_grid]  # Copy the input grid to the output grid

    # Transform the block of 9s
    for i in range(4, 8):
        output_grid[i][3] = 0  # Replace the middle column of 9s with 0s

    # Transform the block of 2s
    for i in range(4, 10):
        output_grid[i][6] = 0  # Replace the middle column of 2s with 0s

    return output_grid

# Test input grid
input_grid = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 9, 9, 9, 2, 2, 2, 2, 0, 0],
    [0, 0, 9, 9, 9, 2, 2, 2, 2, 0, 0],
    [0, 0, 9, 9, 9, 2, 2, 2, 2, 0, 0],
    [0, 0, 9, 9, 9, 2, 2, 2, 2, 0, 0],
    [0, 0, 0, 0, 0, 2, 2, 2, 2, 0, 0],
    [0, 0, 0, 0, 0, 2, 2, 2, 2, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]

output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))