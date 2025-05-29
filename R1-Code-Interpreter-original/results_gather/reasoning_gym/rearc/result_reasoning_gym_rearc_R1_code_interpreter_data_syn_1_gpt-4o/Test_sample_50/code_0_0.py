def transform_grid(input_grid):
    # Find the row with the smallest non-zero number
    min_row = min(input_grid, key=lambda row: min(filter(lambda x: x != 0, row)))
    # Create the output grid
    output_grid = [min_row] + [row for row in input_grid if row != min_row]
    return output_grid

# Test input grid
input_grid = [
    [1, 8, 1],
    [1, 1, 1],
    [1, 1, 1]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))