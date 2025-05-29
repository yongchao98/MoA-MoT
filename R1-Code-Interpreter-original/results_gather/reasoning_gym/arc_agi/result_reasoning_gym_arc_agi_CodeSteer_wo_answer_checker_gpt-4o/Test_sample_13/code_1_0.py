def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Copy the input grid

    for i in range(1, rows - 1):
        for j in range(1, cols - 1):
            # Check if the current '4' is surrounded by '8's in a cross pattern
            if (input_grid[i][j] == 4 and
                input_grid[i-1][j] == 8 and
                input_grid[i+1][j] == 8 and
                input_grid[i][j-1] == 8 and
                input_grid[i][j+1] == 8):
                output_grid[i][j] = 5

    return output_grid

# Test input grid
input_grid = [
    [4, 4, 4, 4, 4, 4, 8, 8, 8, 8, 8, 4, 4, 4, 4, 4],
    [8, 8, 8, 8, 4, 4, 8, 4, 4, 8, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 8, 4, 4, 4, 4, 4, 8, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 8, 8, 8, 8, 8, 8, 8, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    [8, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    [8, 8, 8, 8, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    [8, 4, 4, 8, 4, 4, 4, 4, 8, 8, 8, 8, 8, 8, 4, 4],
    [8, 4, 4, 8, 4, 4, 4, 4, 4, 4, 8, 4, 4, 8, 4, 4],
    [8, 4, 4, 8, 4, 4, 4, 4, 4, 4, 8, 4, 4, 8, 4, 4],
    [8, 4, 4, 8, 4, 4, 4, 4, 4, 4, 8, 8, 8, 8, 4, 4],
    [8, 8, 8, 8, 4, 4, 8, 8, 8, 4, 4, 4, 4, 8, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 8, 4, 4, 8, 4, 4, 4, 4, 4, 4, 8],
    [4, 4, 4, 4, 4, 8, 8, 8, 8, 4, 4, 4, 4, 4, 4, 8]
]

output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))