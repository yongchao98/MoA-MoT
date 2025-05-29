def transform_grid(input_grid):
    # Copy the input grid to the output grid
    output_grid = [row[:] for row in input_grid]
    
    # Define the transformation rule
    for i in range(len(input_grid)):
        for j in range(len(input_grid[i])):
            if input_grid[i][j] == 0:
                # Surround with 1 and replace 0 with 2
                if i > 0:
                    output_grid[i-1][j] = 1
                if i < len(input_grid) - 1:
                    output_grid[i+1][j] = 1
                if j > 0:
                    output_grid[i][j-1] = 1
                if j < len(input_grid[i]) - 1:
                    output_grid[i][j+1] = 1
                output_grid[i][j] = 2
    
    return output_grid

# Test input grid
input_grid = [
    [7, 7, 7, 7, 7, 7, 3],
    [7, 7, 7, 7, 7, 7, 3],
    [7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7],
    [3, 7, 7, 7, 0, 7, 7],
    [3, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))