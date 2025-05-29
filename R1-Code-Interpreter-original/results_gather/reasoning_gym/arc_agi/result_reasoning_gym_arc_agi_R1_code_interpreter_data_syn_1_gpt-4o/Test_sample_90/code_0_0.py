def transform_grid(input_grid):
    # Remove the border
    inner_grid = [row[1:-1] for row in input_grid[1:-1]]
    return inner_grid

# Test input grid
input_grid = [
    [9, 9, 9, 9, 9, 9],
    [9, 7, 5, 7, 5, 7],
    [9, 5, 7, 5, 7, 5],
    [9, 7, 5, 7, 5, 7],
    [9, 5, 7, 5, 7, 5],
    [9, 7, 5, 7, 5, 7]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))