def transform_grid(input_grid):
    # Copy the input grid to the output grid
    output_grid = [row[:] for row in input_grid]
    
    # Define the transformation rules based on observations
    for i in range(len(input_grid)):
        for j in range(len(input_grid[i])):
            if input_grid[i][j] == 1:
                # Replace '1' with '4' based on the pattern observed
                output_grid[i][j] = 4
            elif input_grid[i][j] == 4:
                # Replace '4' with '2' based on the pattern observed
                output_grid[i][j] = 2
    
    return output_grid

# Test input grid
input_grid = [
    [7, 7, 7, 7, 7, 7, 1, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 1, 7, 7, 7, 7, 7],
    [7, 4, 4, 4, 4, 7, 7, 7, 7, 7, 7, 7],
    [7, 4, 7, 2, 4, 7, 7, 7, 7, 7, 7, 7],
    [7, 4, 7, 2, 4, 7, 7, 7, 7, 7, 7, 7],
    [7, 4, 4, 4, 4, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 1, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 1, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 1, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 1, 7, 7, 7]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))