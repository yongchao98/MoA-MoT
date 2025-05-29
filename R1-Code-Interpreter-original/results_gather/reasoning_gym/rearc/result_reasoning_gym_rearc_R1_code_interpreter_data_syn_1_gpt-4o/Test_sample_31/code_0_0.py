def transform_grid(input_grid):
    # Copy the input grid to the output grid
    output_grid = [row[:] for row in input_grid]
    
    # Identify the 3x3 subgrid to transform
    # Based on the examples, it seems to be around the column of 9s
    for i in range(1, len(input_grid) - 1):
        for j in range(1, len(input_grid[0]) - 1):
            if input_grid[i][j] == 9:
                # Transform the 3x3 subgrid
                for x in range(i-1, i+2):
                    for y in range(j-1, j+2):
                        if input_grid[x][y] != 9:
                            output_grid[x][y] = 4
                # Center element becomes 4 if it is 9
                output_grid[i][j] = 4
                return output_grid

# Test input grid
input_grid = [
    [7, 7, 7, 9, 7, 7, 7, 7],
    [7, 7, 7, 9, 7, 7, 7, 7],
    [1, 1, 1, 9, 1, 1, 1, 1],
    [7, 7, 7, 9, 7, 7, 7, 7],
    [7, 7, 7, 9, 7, 7, 7, 7],
    [7, 7, 7, 9, 7, 7, 7, 7],
    [7, 7, 7, 9, 7, 7, 7, 7]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))