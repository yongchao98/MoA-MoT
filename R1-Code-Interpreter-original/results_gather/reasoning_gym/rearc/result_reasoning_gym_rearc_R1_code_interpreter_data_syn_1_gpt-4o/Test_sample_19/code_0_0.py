def transform_grid(input_grid):
    output_grid = [row[:] for row in input_grid]  # Copy the input grid to output grid
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    for i in range(1, rows - 1):
        for j in range(1, cols - 1):
            # Check horizontally
            if input_grid[i][j-1] == input_grid[i][j+1] and input_grid[i][j-1] != input_grid[i][j]:
                output_grid[i][j] = 2
            # Check vertically
            if input_grid[i-1][j] == input_grid[i+1][j] and input_grid[i-1][j] != input_grid[i][j]:
                output_grid[i][j] = 2
    
    return output_grid

# Test input grid
input_grid = [
    [7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7],
    [8, 7, 8, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7]
]

output_grid = transform_grid(input_grid)

for row in output_grid:
    print(' '.join(map(str, row)))