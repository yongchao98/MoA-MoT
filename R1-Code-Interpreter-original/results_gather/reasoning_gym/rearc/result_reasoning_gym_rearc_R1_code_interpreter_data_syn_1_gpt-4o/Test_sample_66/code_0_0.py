def transform_grid(input_grid):
    output_grid = [row[:] for row in input_grid]  # Copy the input grid to the output grid

    # Define the transformation rules based on the examples
    for i in range(len(input_grid)):
        for j in range(len(input_grid[i])):
            if input_grid[i][j] == 7:
                # Check the context of the number 7
                if (i > 0 and input_grid[i-1][j] == 7) or (j > 0 and input_grid[i][j-1] == 7):
                    output_grid[i][j] = 8
                else:
                    output_grid[i][j] = 2

    return output_grid

# Test input grid
input_grid = [
    [1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 7, 7],
    [1, 1, 1, 1, 1, 1, 1, 7, 7],
    [1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 7, 7, 7, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1]
]

output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))