def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Copy the input grid

    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] not in (0, 9, 6):
                # Check horizontal sequence
                if (j > 0 and input_grid[i][j - 1] not in (0, 9, 6)) or (j < cols - 1 and input_grid[i][j + 1] not in (0, 9, 6)):
                    output_grid[i][j] = 5
                # Check vertical sequence
                if (i > 0 and input_grid[i - 1][j] not in (0, 9, 6)) or (i < rows - 1 and input_grid[i + 1][j] not in (0, 9, 6)):
                    output_grid[i][j] = 5

    return output_grid

# Test input grid
input_grid = [
    [3, 9, 3, 3, 0, 3, 2, 4, 3, 3, 6, 3, 3, 3, 6],
    [3, 3, 3, 5, 3, 3, 7, 3, 3, 2, 6, 3, 3, 0, 0],
    [3, 6, 3, 3, 3, 6, 0, 3, 2, 7, 3, 2, 3, 3, 3],
    [7, 3, 7, 3, 0, 0, 3, 5, 0, 3, 3, 3, 3, 3, 3],
    [3, 4, 2, 9, 3, 0, 3, 3, 2, 7, 3, 5, 2, 3, 3],
    [3, 7, 3, 3, 5, 3, 3, 3, 0, 2, 4, 3, 3, 3, 3],
    [0, 3, 0, 3, 3, 3, 3, 4, 7, 3, 3, 2, 3, 0, 2],
    [3, 3, 9, 9, 3, 6, 7, 3, 0, 7, 3, 4, 0, 4, 3],
    [7, 7, 7, 0, 3, 2, 7, 3, 5, 3, 9, 7, 3, 5, 3],
    [3, 3, 3, 3, 9, 3, 6, 3, 0, 2, 3, 0, 4, 3, 3],
    [3, 7, 5, 3, 3, 3, 6, 3, 0, 3, 3, 7, 3, 4, 9],
    [3, 3, 3, 3, 3, 2, 2, 7, 3, 3, 4, 2, 3, 3, 9],
    [3, 3, 9, 3, 2, 3, 3, 3, 6, 2, 2, 3, 4, 3, 3],
    [3, 3, 7, 7, 0, 2, 3, 2, 3, 3, 6, 2, 3, 5, 3],
    [4, 2, 9, 2, 6, 7, 7, 3, 7, 6, 2, 7, 3, 4, 0]
]

output_grid = transform_grid(input_grid)

for row in output_grid:
    print(' '.join(map(str, row)))