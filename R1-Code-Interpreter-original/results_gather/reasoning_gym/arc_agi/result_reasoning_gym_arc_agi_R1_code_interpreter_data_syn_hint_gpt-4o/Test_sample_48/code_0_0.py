def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Copy the input grid

    for i in range(1, rows - 1):
        for j in range(1, cols - 1):
            if input_grid[i][j] == 1:
                # Check if the current '1' is surrounded by other '1's
                if (input_grid[i-1][j] == 1 and input_grid[i+1][j] == 1) or \
                   (input_grid[i][j-1] == 1 and input_grid[i][j+1] == 1):
                    output_grid[i][j] = 5

    return output_grid

# Test input grid
input_grid = [
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 1, 1],
    [6, 6, 6, 1, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 1, 6, 6, 6, 1, 1, 1, 1, 1, 6, 6, 6],
    [6, 6, 6, 1, 6, 6, 6, 1, 6, 6, 1, 6, 6, 6, 6],
    [6, 6, 6, 1, 6, 6, 6, 1, 6, 6, 1, 6, 6, 6, 6],
    [1, 6, 6, 1, 6, 6, 6, 1, 1, 1, 1, 6, 6, 6, 6],
    [1, 1, 1, 1, 6, 6, 6, 1, 6, 6, 6, 6, 6, 6, 6],
    [1, 6, 6, 1, 6, 6, 6, 1, 6, 6, 6, 1, 1, 1, 1],
    [1, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 1, 6, 6, 1],
    [1, 1, 6, 6, 6, 6, 6, 6, 6, 6, 6, 1, 6, 6, 1],
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 1, 1],
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 1, 1, 6, 6, 6, 1, 1, 1, 1, 1, 1, 6, 6, 6],
    [6, 1, 6, 6, 6, 6, 1, 6, 6, 6, 6, 1, 6, 6, 6],
    [6, 1, 6, 6, 6, 6, 1, 6, 6, 6, 6, 1, 6, 6, 6],
    [6, 1, 6, 6, 6, 1, 1, 1, 1, 1, 1, 1, 6, 6, 6]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))