def transform_grid(input_grid):
    # Copy the input grid to output grid
    output_grid = [row[:] for row in input_grid]
    
    # Apply the transformation rule
    # Replace the first row with the last row
    output_grid[0] = input_grid[-1]
    
    return output_grid

# Test input
input_grid = [
    [7, 7, 0, 7, 7],
    [2, 7, 0, 7, 7]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))