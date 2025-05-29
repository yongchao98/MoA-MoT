def transform_grid(input_grid):
    # Repeat the input grid twice vertically
    output_grid = input_grid + input_grid
    return output_grid

# Test input grid
input_grid = [
    [8, 8, 8],
    [8, 8, 8]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))