def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Copy the input grid

    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] == 0:
                # Check adjacent cells
                if i > 0 and input_grid[i-1][j] != 0:  # Up
                    output_grid[i-1][j] = 5
                if i < rows - 1 and input_grid[i+1][j] != 0:  # Down
                    output_grid[i+1][j] = 5
                if j > 0 and input_grid[i][j-1] != 0:  # Left
                    output_grid[i][j-1] = 5
                if j < cols - 1 and input_grid[i][j+1] != 0:  # Right
                    output_grid[i][j+1] = 5

    return output_grid

# Test input grid
input_grid = [
    [1, 7, 1, 1, 1, 1, 1, 1, 7, 7, 1, 1, 7, 7, 7, 1, 1, 1],
    [1, 1, 1, 1, 1, 7, 1, 1, 7, 1, 1, 1, 7, 0, 0, 7, 0, 7],
    [1, 7, 1, 1, 7, 1, 1, 7, 7, 7, 1, 1, 0, 7, 7, 0, 0, 1],
    [1, 7, 1, 7, 7, 7, 7, 7, 7, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [7, 7, 7, 7, 1, 1, 7, 1, 7, 7, 1, 7, 1, 1, 1, 7, 7, 1],
    [1, 1, 1, 1, 7, 1, 1, 7, 0, 0, 7, 0, 7, 1, 7, 7, 7, 7],
    [1, 1, 1, 1, 7, 7, 1, 7, 0, 0, 0, 0, 7, 1, 1, 1, 1, 1],
    [7, 1, 7, 7, 7, 1, 7, 7, 1, 7, 1, 1, 7, 1, 7, 7, 1, 7],
    [1, 1, 1, 1, 7, 1, 7, 7, 7, 7, 7, 1, 1, 1, 7, 1, 1, 1],
    [7, 1, 1, 1, 1, 7, 7, 1, 1, 7, 7, 1, 7, 1, 7, 1, 7, 7],
    [7, 7, 7, 7, 1, 1, 7, 1, 1, 7, 0, 7, 1, 7, 1, 7, 7, 1],
    [7, 1, 7, 7, 7, 1, 1, 1, 7, 7, 7, 0, 7, 1, 1, 7, 7, 1],
    [1, 7, 7, 7, 1, 1, 1, 7, 1, 1, 7, 7, 7, 1, 7, 7, 1, 7],
    [1, 7, 0, 0, 7, 0, 0, 7, 1, 1, 1, 1, 7, 7, 7, 7, 7, 7],
    [7, 0, 7, 7, 0, 7, 0, 1, 1, 7, 1, 7, 7, 7, 7, 7, 7, 1],
    [1, 7, 0, 0, 7, 0, 0, 1, 7, 1, 7, 7, 1, 1, 7, 1, 7, 1],
    [1, 7, 7, 7, 1, 7, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))