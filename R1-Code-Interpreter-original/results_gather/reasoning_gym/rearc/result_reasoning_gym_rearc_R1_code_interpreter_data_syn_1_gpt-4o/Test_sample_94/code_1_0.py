def transform_grid(input_grid):
    # Create a copy of the input grid to modify
    output_grid = [row[:] for row in input_grid]
    
    # Replace the first row with the second row
    if len(input_grid) > 1:
        output_grid[0] = input_grid[1][:]
        
        # Ensure the third element of the first row remains unchanged
        if len(input_grid[0]) > 2:
            output_grid[0][2] = input_grid[0][2]
    
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