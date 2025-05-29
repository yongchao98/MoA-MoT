def transform_grid(input_grid):
    # Copy the input grid to the output grid
    output_grid = [row[:] for row in input_grid]
    
    # Apply the rule to the first row
    if len(input_grid[0]) > 1:
        output_grid[0][-2] = input_grid[0][3]  # Place the number at the second last position
    
    return output_grid

# Test input grid
input_grid = [
    [9, 9, 8, 2, 9, 9, 9],
    [6, 9, 9, 2, 9, 9, 9],
    [9, 9, 9, 2, 9, 9, 9]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))