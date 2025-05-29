def transform_grid(input_grid):
    output_grid = [row[:] for row in input_grid]  # Copy the input grid to the output grid

    rows = len(input_grid)
    cols = len(input_grid[0])

    # Check for 2x2 blocks of '1's at the edges and transform them to '0's
    for i in range(rows - 1):
        for j in range(cols - 1):
            if input_grid[i][j] == 1 and input_grid[i][j + 1] == 1 and input_grid[i + 1][j] == 1 and input_grid[i + 1][j + 1] == 1:
                # Check if the block is at the top, bottom, left, or right edge
                if i == 0 or j == 0 or i == rows - 2 or j == cols - 2:
                    output_grid[i][j] = 0
                    output_grid[i][j + 1] = 0
                    output_grid[i + 1][j] = 0
                    output_grid[i + 1][j + 1] = 0

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

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))