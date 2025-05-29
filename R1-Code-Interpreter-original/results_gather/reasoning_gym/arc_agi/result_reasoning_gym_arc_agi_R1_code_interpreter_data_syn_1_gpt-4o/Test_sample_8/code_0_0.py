def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Copy the input grid

    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] == 8:
                # Check and replace adjacent cells
                if i > 0 and input_grid[i-1][j] != 8:
                    output_grid[i-1][j] = 0
                if i < rows - 1 and input_grid[i+1][j] != 8:
                    output_grid[i+1][j] = 0
                if j > 0 and input_grid[i][j-1] != 8:
                    output_grid[i][j-1] = 0
                if j < cols - 1 and input_grid[i][j+1] != 8:
                    output_grid[i][j+1] = 0

    return output_grid

# Test input grid
input_grid = [
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 8, 5, 8, 8, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 8, 5, 5, 8, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 8, 5, 5, 8, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 8, 8, 8, 8, 5, 5, 5, 5],
    [5, 8, 8, 8, 8, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 8, 5, 5, 8, 5, 5, 5, 8, 8, 8, 8, 5, 5, 5],
    [5, 8, 5, 5, 8, 5, 5, 5, 8, 5, 5, 8, 5, 5, 5],
    [5, 8, 5, 8, 8, 5, 5, 5, 8, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 8, 8, 8, 8, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))