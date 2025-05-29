def transform_grid(input_grid):
    # Determine the size of the grid
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Create an output grid filled with 8s
    output_grid = [[8 for _ in range(cols)] for _ in range(rows)]
    
    # Preserve the last few rows and columns
    for i in range(rows - 5, rows):
        for j in range(cols - 5, cols):
            output_grid[i][j] = input_grid[i][j]
    
    return output_grid

# Test input grid
input_grid = [
    [8, 8, 8, 8, 8, 8, 8, 8, 6, 6, 6, 8],
    [8, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 8],
    [8, 6, 9, 6, 6, 6, 6, 8, 8, 8, 8, 8],
    [8, 6, 6, 6, 6, 9, 6, 8, 6, 6, 6, 6],
    [8, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6],
    [8, 8, 8, 8, 8, 8, 8, 8, 6, 6, 9, 6],
    [8, 8, 8, 8, 6, 6, 6, 8, 6, 6, 6, 6],
    [8, 8, 8, 8, 6, 6, 9, 8, 8, 8, 8, 8],
    [8, 8, 8, 8, 6, 6, 6, 8, 8, 8, 8, 8],
    [6, 6, 6, 8, 9, 6, 6, 8, 8, 8, 8, 8],
    [9, 6, 6, 8, 8, 8, 8, 8, 6, 9, 6, 9],
    [6, 6, 6, 8, 8, 8, 8, 8, 6, 6, 6, 6],
    [8, 8, 8, 8, 8, 8, 8, 8, 6, 9, 6, 6]
]

# Transform the input grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))