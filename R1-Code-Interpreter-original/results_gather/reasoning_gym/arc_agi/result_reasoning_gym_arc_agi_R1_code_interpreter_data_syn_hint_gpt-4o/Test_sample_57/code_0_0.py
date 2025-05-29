def transform_grid(input_grid):
    n = len(input_grid)
    m = len(input_grid[0])
    
    # Create a new grid for the output
    output_grid = [row[:] for row in input_grid]
    
    # Apply the diagonal increment rule
    for i in range(n - 1):
        output_grid[i][i] = 9
    
    # Apply the last row transformation
    for j in range(m - 1):
        output_grid[n - 1][j] = 1
    
    return output_grid

# Test input grid
input_grid = [
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 3],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 3],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 3],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 3],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 3],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 3],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 3],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 3],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 3],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 3]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(row)