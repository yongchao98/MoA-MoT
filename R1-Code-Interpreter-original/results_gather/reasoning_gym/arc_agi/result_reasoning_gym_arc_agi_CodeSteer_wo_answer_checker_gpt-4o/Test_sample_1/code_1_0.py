def transform_grid(input_grid):
    # Initialize an empty 8x8 grid
    output_grid = [[0] * 8 for _ in range(8)]
    
    # Fill the output grid based on observed patterns
    # This pattern is derived from analyzing the examples
    for i in range(4):
        for j in range(4):
            output_grid[i][j] = input_grid[i][j]
            output_grid[i][j+4] = input_grid[j][3-i]
            output_grid[i+4][j] = input_grid[3-j][i]
            output_grid[i+4][j+4] = input_grid[3-i][3-j]
    
    return output_grid

# Test input grid
input_grid = [
    [4, 6, 4, 2],
    [2, 5, 2, 4],
    [6, 5, 6, 6],
    [6, 5, 6, 6]
]

# Get the output grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))