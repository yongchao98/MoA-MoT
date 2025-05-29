def transform_grid(input_grid):
    # Remove rows with unique patterns (in this case, rows containing 9s)
    output_grid = [row for row in input_grid if 9 not in row]
    return output_grid

# Test input grid
input_grid = [
    [0, 6, 0, 6, 0, 6, 0, 6],
    [6, 0, 6, 0, 6, 0, 6, 0],
    [0, 6, 0, 6, 0, 6, 0, 6],
    [6, 0, 6, 0, 6, 9, 9, 0],
    [0, 6, 0, 6, 0, 9, 9, 6],
    [6, 0, 6, 0, 6, 0, 6, 0],
    [0, 6, 0, 6, 0, 6, 0, 6],
    [6, 0, 6, 0, 6, 0, 6, 0]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))