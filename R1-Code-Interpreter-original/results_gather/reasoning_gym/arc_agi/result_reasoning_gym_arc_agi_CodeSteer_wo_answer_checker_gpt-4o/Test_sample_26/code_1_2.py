def transform_grid(input_grid):
    output_grid = [row[:] for row in input_grid]  # Copy the input grid to the output grid

    rows = len(input_grid)
    cols = len(input_grid[0])

    # Check for 2x2 blocks of '1's and transform them to '0's if they are at or adjacent to the edge
    for i in range(rows - 1):
        for j in range(cols - 1):
            if input_grid[i][j] == 1 and input_grid[i][j + 1] == 1 and input_grid[i + 1][j] == 1 and input_grid[i + 1][j + 1] == 1:
                # Check if the block is at the edge or adjacent to the edge
                if i <= 1 or j <= 1 or i >= rows - 3 or j >= cols - 3:
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