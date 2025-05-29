def transform_grid(input_grid):
    # Copy the input grid to the output grid
    output_grid = [row[:] for row in input_grid]
    
    # Get the number of rows and columns
    num_rows = len(input_grid)
    num_cols = len(input_grid[0])
    
    # Transform the last column
    for i in range(num_rows):
        if input_grid[i][num_cols - 1] == 0:
            output_grid[i][num_cols - 1] = 1
    
    return output_grid

# Test input grid
input_grid = [
    [2, 7, 8, 0, 0, 0, 0],
    [9, 0, 8, 8, 0, 0, 2],
    [0, 0, 3, 0, 0, 7, 0],
    [0, 0, 0, 0, 0, 0, 4],
    [2, 4, 0, 0, 0, 1, 0],
    [3, 0, 0, 6, 2, 0, 0],
    [6, 6, 7, 0, 4, 4, 0]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))