def transform_grid(input_grid):
    output_grid = [row[:] for row in input_grid]  # Copy the input grid to the output grid

    # Transform rows
    for i in range(len(input_grid)):
        for j in range(len(input_grid[i]) - 1):
            if input_grid[i][j] == 1 and input_grid[i][j + 1] == 1:
                output_grid[i][j] = 0
                output_grid[i][j + 1] = 0
                break  # Only transform the first pair

    # Transform columns
    for j in range(len(input_grid[0])):
        for i in range(len(input_grid) - 1):
            if input_grid[i][j] == 1 and input_grid[i + 1][j] == 1:
                output_grid[i][j] = 0
                output_grid[i + 1][j] = 0
                break  # Only transform the first pair

    return output_grid

# Test input grid
input_grid = [
    [8, 8, 8, 8, 8, 1, 1],
    [8, 8, 8, 8, 8, 1, 1],
    [8, 8, 1, 1, 8, 8, 8],
    [8, 8, 1, 1, 8, 8, 8],
    [8, 8, 8, 8, 1, 1, 8],
    [8, 8, 8, 8, 1, 1, 8],
    [8, 8, 8, 8, 8, 8, 8],
    [8, 1, 1, 8, 8, 8, 8],
    [8, 1, 1, 8, 8, 1, 1],
    [8, 8, 8, 8, 8, 1, 1],
    [8, 8, 8, 8, 8, 8, 8],
    [8, 8, 8, 8, 1, 1, 8],
    [8, 8, 8, 8, 1, 1, 8],
    [8, 1, 1, 8, 8, 8, 8],
    [8, 1, 1, 8, 8, 8, 8]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))