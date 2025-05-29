def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Start with a copy of the input grid

    def fill_with_threes(r, c):
        if 0 <= r < rows and 0 <= c < cols and output_grid[r][c] == 7:
            output_grid[r][c] = 3
            fill_with_threes(r + 1, c)
            fill_with_threes(r - 1, c)
            fill_with_threes(r, c + 1)
            fill_with_threes(r, c - 1)

    # Find all '3's in the input grid and start filling from there
    for r in range(rows):
        for c in range(cols):
            if input_grid[r][c] == 3:
                fill_with_threes(r, c)

    # Check for regions surrounded by '1's and fill them with '3's
    for r in range(rows):
        for c in range(cols):
            if input_grid[r][c] == 1:
                if r > 0 and input_grid[r - 1][c] == 7:
                    fill_with_threes(r - 1, c)
                if r < rows - 1 and input_grid[r + 1][c] == 7:
                    fill_with_threes(r + 1, c)
                if c > 0 and input_grid[r][c - 1] == 7:
                    fill_with_threes(r, c - 1)
                if c < cols - 1 and input_grid[r][c + 1] == 7:
                    fill_with_threes(r, c + 1)

    return output_grid

# Test input grid
input_grid = [
    [7, 7, 7, 1, 7, 7, 7, 1, 7, 7, 7, 1, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 1, 1, 1, 7, 7],
    [7, 7, 7, 1, 1, 7, 1, 7, 7, 7, 7, 1, 7],
    [1, 7, 7, 7, 1, 7, 1, 7, 7, 1, 7, 7, 3],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 1, 7],
    [1, 1, 7, 7, 1, 7, 1, 1, 7, 7, 7, 7, 7],
    [7, 7, 1, 1, 7, 7, 7, 1, 7, 7, 7, 1, 7],
    [7, 1, 7, 1, 7, 7, 1, 7, 7, 1, 7, 7, 3],
    [7, 7, 7, 7, 7, 1, 7, 1, 7, 7, 7, 7, 7],
    [7, 1, 7, 7, 7, 7, 7, 7, 1, 1, 1, 7, 3],
    [7, 1, 7, 7, 7, 1, 7, 7, 1, 7, 1, 7, 1],
    [1, 1, 7, 7, 7, 7, 7, 7, 7, 1, 7, 7, 1],
    [1, 7, 7, 1, 1, 7, 7, 7, 7, 7, 7, 7, 7]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))