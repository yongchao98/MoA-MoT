def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Copy the input grid

    for r in range(rows):
        for c in range(cols):
            if input_grid[r][c] == 4:
                # Check left for adding 3
                if c > 0 and input_grid[r][c-1] == 0:
                    output_grid[r][c-1] = 3
                # Check above for adding 1
                if r > 0 and input_grid[r-1][c] == 0:
                    output_grid[r-1][c] = 1
                # Check below for adding 1
                if r < rows - 1 and input_grid[r+1][c] == 0:
                    output_grid[r+1][c] = 1

    return output_grid

# Test input grid
input_grid = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0],
    [0, 0, 4, 4, 0, 0, 0, 0, 0, 0, 6, 4, 4, 0, 0],
    [0, 0, 4, 0, 4, 0, 0, 0, 0, 1, 4, 0, 4, 0, 0],
    [0, 0, 6, 4, 4, 0, 0, 0, 0, 3, 4, 4, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 4, 4, 6, 0, 0, 0, 0, 0, 4, 4, 0],
    [0, 0, 0, 0, 4, 0, 4, 0, 0, 0, 0, 4, 0, 4, 0],
    [0, 0, 0, 0, 0, 4, 4, 0, 0, 0, 0, 4, 4, 6, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))