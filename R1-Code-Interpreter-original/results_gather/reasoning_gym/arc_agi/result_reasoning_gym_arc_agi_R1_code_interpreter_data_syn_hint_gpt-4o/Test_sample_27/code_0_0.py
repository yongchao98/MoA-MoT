def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Copy the input grid

    for r in range(rows):
        for c in range(cols):
            if input_grid[r][c] == 9:
                # Add 5 to the left if possible
                if c > 0 and input_grid[r][c-1] != 9:
                    output_grid[r][c-1] = 5
                # Add 4 to the right if possible
                if c < cols - 1 and input_grid[r][c+1] != 9:
                    output_grid[r][c+1] = 4
                # Add 1 below if possible
                if r < rows - 1 and input_grid[r+1][c] != 9:
                    output_grid[r+1][c] = 1
                # Add 8 above if possible
                if r > 0 and input_grid[r-1][c] != 9:
                    output_grid[r-1][c] = 8

    return output_grid

# Test input grid
input_grid = [
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    [3, 9, 9, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    [3, 9, 9, 3, 3, 3, 3, 3, 3, 3, 9, 9, 3, 3, 3],
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 9, 9, 3, 3, 3],
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    [3, 3, 3, 3, 9, 9, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    [3, 3, 3, 3, 9, 9, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 9, 9, 3, 3],
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 9, 9, 3, 3],
    [3, 9, 9, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    [3, 9, 9, 3, 3, 3, 9, 9, 3, 3, 3, 3, 3, 3, 3],
    [3, 3, 3, 3, 3, 3, 9, 9, 3, 3, 3, 3, 3, 3, 3],
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
]

output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))