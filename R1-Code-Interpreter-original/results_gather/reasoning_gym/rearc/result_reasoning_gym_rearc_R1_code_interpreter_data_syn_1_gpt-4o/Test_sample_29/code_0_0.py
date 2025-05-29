def transform_grid(input_grid):
    # Copy the input grid to the output grid
    output_grid = [row[:] for row in input_grid]
    
    # Define the transformation rule based on the examples
    # Here, we replace '8' with alternating '7' and '1' in the lower half of the grid
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Start replacing from the middle row downwards
    for i in range(rows // 2, rows):
        for j in range(cols):
            if input_grid[i][j] == 8:
                # Alternate between 7 and 1
                output_grid[i][j] = 7 if (i + j) % 2 == 0 else 1
    
    return output_grid

# Test input grid
input_grid = [
    [1, 1, 1, 1, 1, 8, 1],
    [8, 8, 1, 2, 7, 8, 1],
    [1, 1, 1, 7, 1, 1, 1],
    [1, 8, 8, 8, 1, 8, 8],
    [1, 1, 1, 1, 1, 1, 1]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))